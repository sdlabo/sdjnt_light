#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include <time.h>
#include <float.h>

#include "fftsg.h"
#include "sdlab.h"

#define DATA_SIZE (2*1024*1024)  // 2^(21) points
#define PORTA (0x4000) // receive port #.
// the size of a buffer for receiving a UDP packet
#define RECV_BUF_SIZE   (2048)
// the size of a length of data in a UDP packet
#define DATA_BURST_SIZE  (512)
// the buffer size more than 2 has not been supported yet
#define NUM_OF_ENV_BUFFER (2)

// a working data for receiving packet and FFT operation
struct env{
  // a pointer to the instance of UDP communication
  int sock;
  int port;
  // receiving data buffers, selecting as double buffer
  double *buf1[NUM_OF_ENV_BUFFER];
  double *buf2[NUM_OF_ENV_BUFFER];
  // a pointer to the active buffer to receive UDP packet.
  double *cur1;
  double *cur2;
  // a pointer to the buffer working FFT operation
  double *fft_work1;
  double *fft_work2;
  // internal data for FFT operation
  int *ip1;
  int *ip2;
  // internal data for FFT operation
  double *w1;
  double *w2;
  // temporal buffer to receive a UDP packet.
  char recv_buf[RECV_BUF_SIZE];
  // identifier for the active buffer;
  int buf_id1;
  int buf_id2;
};

struct fft_arg
{
  double *fft_work;
  int *ip;
  double *w;
};

extern char param_basename[1024];
extern int param_duration;
extern int param_accumulation_time;
extern bool param_is_output_rawfile;

bool g_is_calc_finished;

// a system envrionment
struct system_env{
  int file_count; // postfix #. of file to stored result data
  // result buffers to save results into HDD in a group, used as double buffer
  double * result[2];
  double * cur; // a pointer to the active buffer to pack result
  double * dump; // a pointer to the buffer to save HDD
  int result_count; // the number of buffering results
};

void *fft_thread(void *param); // FFT wrapper function.

int udp_init(in_addr_t in_addr, int port)
{
  int sock;
  struct sockaddr_in addr;

  sock = socket(AF_INET, SOCK_DGRAM, 0);
  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  addr.sin_addr.s_addr = in_addr;

  int n = 16 * 1024 * 1024;
  if(setsockopt(sock, SOL_SOCKET, SO_RCVBUF, &n, sizeof(n)) == -1) {
    perror("setsockopt");
    exit(1);
  }

  bind(sock, (struct sockaddr*)&addr, sizeof(struct sockaddr_in));

  return sock;
}


static void * malloc_or_die(size_t sz, const char *msg)
{
  void* ptr = malloc(sz);
  if(ptr == NULL){
    fprintf(stderr, "%s\n", msg);
    exit(0);
  }
  return ptr;
}

static void sys_init(struct system_env *sys){
  sys->result[0] = (double*)malloc_or_die(
    sizeof(double) * DATA_SIZE / 2 * 4,
    "cannot allocate save result buffer (0)");
  sys->result[1] = (double*)malloc_or_die(
    sizeof(double) * DATA_SIZE /2 * 4,
    "cannot allocate save result buffer (1)");

  sys->result_count = 0;
  sys->file_count = 0;
  sys->cur = sys->result[0];
}

void sys_free(struct system_env *sys){
  free(sys->result[0]);
  free(sys->result[1]);
}

void env_init(struct env *e, int port){
  // required memories are allocated
  for(int i = 0; i < NUM_OF_ENV_BUFFER; i++){
    e->buf1[i] = (double*)malloc_or_die(
      sizeof(double) * DATA_SIZE * 2,
      "cannot allocate heap memory for I/O data");
  }

  for(int i = 0; i < NUM_OF_ENV_BUFFER; i++){
    e->buf2[i] = (double*)malloc_or_die(
      sizeof(double) * DATA_SIZE * 2,
      "cannot allocate heap memory for I/O data");
  }

  e->ip1 = (int*)malloc_or_die(
    sizeof(int) * (2 + sqrt(DATA_SIZE*2)),
    "cannot allocate heap memory for bit reversal");

  e->ip2 = (int*)malloc_or_die(
    sizeof(int) * (2 + sqrt(DATA_SIZE*2)),
    "cannot allocate heap memory for bit reversal");

  e->w1  = (double*)malloc_or_die(
    sizeof(double) * DATA_SIZE,
    "cannot allocate heap memory for cos/sin table");

  e->w2  = (double*)malloc_or_die(
    sizeof(double) * DATA_SIZE,
    "cannot allocate heap memory for cos/sin table");

  // poiters are initialized.
  e->buf_id1   = 0;
  e->buf_id2   = 0;

  e->cur1      = e->buf1[e->buf_id1];
  e->cur2      = e->buf2[e->buf_id2];

  e->fft_work1 = e->buf1[e->buf_id1];
  e->fft_work2 = e->buf2[e->buf_id2];

  // to make the constant table for FFT before actual operation.
  e->ip1[0] = 0;
  e->ip2[0] = 0;

  struct fft_arg fe1, fe2;

  fe1.fft_work = e->fft_work1;
  fe1.ip = e->ip1;
  fe1.w = e->w1;

  fe2.fft_work = e->fft_work2;
  fe2.ip = e->ip2;
  fe2.w = e->w2;

  fft_thread(&fe1);
  fft_thread(&fe2);

  e->sock = udp_init(INADDR_ANY, port);
  e->port = port;
}

void env_free(struct env *e){
  free(e->buf1[0]);
  free(e->buf1[1]);
  free(e->ip1);
  free(e->w1);

  free(e->buf2[0]);
  free(e->buf2[1]);
  free(e->ip2);
  free(e->w2);
}

void *fft_thread(void *param)
{
  struct fft_arg *fft = (struct fft_arg*) param;

  cdft(DATA_SIZE * 2, -1, fft->fft_work, fft->ip, fft->w);

  return NULL;
}



struct calc_arg{
  struct system_env *sys;
  struct env *e; // pointer for F(a)
};

/**
 * save results of F(a) and F(b) into the internal buffer.
 * data format in the internal buffer is as fowllows;
 *  - cross-Re
 *  - cross-Im
 *  - F(a)-Re
 *  - F(b)-Re
 */
void save_result(struct calc_arg *arg){
  double a, b, c, d;
  double x_re[DATA_SIZE / 2];
  double x_im[DATA_SIZE / 2];
  double p_a[DATA_SIZE / 2];
  double p_b[DATA_SIZE / 2];

  for(int i = 0; i < DATA_SIZE / 2; i++){
    a = arg->e->fft_work1[2 * i];   // F(a)-Re
    b = arg->e->fft_work1[2 * i + 1]; // F(a)-Im
    c = arg->e->fft_work2[2 * i];   // F(b)-Re
    d = arg->e->fft_work2[2 * i + 1]; // F(b)-Im
    x_re[i] = a * c + b * d;
    x_im[i] = b * c - a * d;
    p_a[i] = a * a + b * b;
    p_b[i] = c * c + d * d;
    arg->sys->cur[DATA_SIZE / 2 * 0 + i] += x_re[i];
    arg->sys->cur[DATA_SIZE / 2 * 1 + i] += x_im[i];
    arg->sys->cur[DATA_SIZE / 2 * 2 + i] += p_a[i];
    arg->sys->cur[DATA_SIZE / 2 * 3 + i] += p_b[i];
  }

  // increment the counter for the number of stored results. 
  arg->sys->result_count++;
}

/**
 * save results stored in the internal buffer into HDD in a group.
 */
//void dump_result_to_file(struct calc_arg *arg){
void* dump_thread(void *param){
  struct calc_arg *arg = (struct calc_arg*) param;
  // make file name and open the file
  char str[256];
  sprintf(str, "%s_%08d_%d",
          param_basename,
          arg->sys->file_count,
          param_accumulation_time);

  printf("opening %s\n", str);
  FILE *fp;
  // the file is created newly at first.
  printf("create: %s\n", str);
  fp = fopen(str, "wb");

  if(fp == NULL){
    fprintf(stderr, "cannot open file: %s", str);
    exit(0);
  }

  // write all data in the internal buffer
  fwrite(arg->sys->dump, sizeof(double), 4 * DATA_SIZE / 2, fp);
  fflush(fp);
  fclose(fp);

  bzero(arg->sys->dump, sizeof(double) * 4 * DATA_SIZE / 2);
  arg->sys->file_count++;

  printf("\t\t\tfile save done\n");

  return NULL;
}

/**
 * FFT operation for received data.
 */
void *calc_thread(void *param)
{
  struct calc_arg *arg = (struct calc_arg*) param;
  pthread_t fft1, fft2;

  struct fft_arg fe1, fe2;

  if(param_is_output_rawfile == true){
    FILE* fp = fopen("/tmp/ch1_work.raw", "w");
    fwrite(arg->e->fft_work1, sizeof(double), 2 * DATA_SIZE, fp);
    fclose(fp);

    fp = fopen("/tmp/ch2_work.raw", "w");
    fwrite(arg->e->fft_work2, sizeof(double), 2 * DATA_SIZE, fp);
    fclose(fp);

    rename("/tmp/ch1_work.raw", "/tmp/sdjnt_light_ch1.raw");
    rename("/tmp/ch2_work.raw", "/tmp/sdjnt_light_ch2.raw");
  }

  fe1.fft_work = arg->e->fft_work1;
  fe1.ip = arg->e->ip1;
  fe1.w = arg->e->w1;

  fe2.fft_work = arg->e->fft_work2;
  fe2.ip = arg->e->ip2;
  fe2.w = arg->e->w2;

  pthread_create(&fft1, NULL, fft_thread, &fe1);
  pthread_create(&fft2, NULL, fft_thread, &fe2);

  pthread_join(fft1, NULL);
  pthread_join(fft2, NULL);

  save_result(arg);

  if(arg->sys->result_count == param_accumulation_time){
    arg->sys->result_count = 0;
    arg->sys->dump = arg->sys->cur; // to point the buffer keeping results.
    arg->sys->cur = (arg->sys->cur == arg->sys->result[0]) ?
      arg->sys->result[1] : arg->sys->result[0];
    printf("\t\tkick save\n");

    printf("saving file\n");
    pthread_t th;
    pthread_create(&th, NULL, dump_thread, arg);
    pthread_detach(th);

  }
  printf("\tcalc done\n");
  g_is_calc_finished = true;

  return NULL;
}

/**
 * Receive UDP packets:
 *  Each UDP packet contains DATA_BURST_SIZE data.
 *  All packets for FFT (DATA_SIZE) is received in this function.
 */
void *recv_thread(void *param){
  struct env *e = (struct env*) param;
  int idx = 0;
  int id = 0, prev_id = 0;
  while(idx < DATA_SIZE * 2){
    recv(e->sock, e->recv_buf, sizeof(int) * RECV_BUF_SIZE, 0);

    int *pi = (int*) e->recv_buf;
    id = ntohl(*pi);
    if(idx == 0 && id != 0){
      // drop packets until receiving a start packet.
      continue;
    }

    if((id - prev_id) > DATA_BURST_SIZE){
      printf("!!! Drop a packet: %08x - %08x = %08x\n",
             id, prev_id, id - prev_id);
    }

    prev_id = id;
    short *data = (short*)(e->recv_buf + sizeof(int));
    for(int i = 0; i < DATA_BURST_SIZE; i++){
      short s = ntohs(data[i]);
      if(i % 2 == 0){
        e->cur1[idx + i + 0] = ((double) s); // Re
        e->cur1[idx + i + 1] = (double) 0; // Im
      }else{
        e->cur2[idx + i - 1 + 0] = ((double) s); // Re
        e->cur2[idx + i - 1 + 1] = (double) 0; // Im
      }
    }

    idx += DATA_BURST_SIZE;
  }


  return NULL;
}

/**
 * swap the working buffer in round-robin manner.
 */
void swap_buffer(struct env *e){
  e->fft_work1 = e->cur1;
  e->buf_id1 = (e->buf_id1 == NUM_OF_ENV_BUFFER-1) ? 0 : e->buf_id1 + 1;
  e->cur1 = e->buf1[e->buf_id1];

  e->fft_work2 = e->cur2;
  e->buf_id2 = (e->buf_id2 == NUM_OF_ENV_BUFFER-1) ? 0 : e->buf_id2 + 1;
  e->cur2 = e->buf2[e->buf_id2];
}


void* sdlab_signal_main()
{
  struct system_env sys;
  struct env e;
  struct calc_arg arg;

  sys_init(&sys);
  env_init(&e, PORTA);

  arg.sys = &sys;
  arg.e = &e;

  pthread_attr_t tattr;
  pthread_attr_init(&tattr);
  pthread_attr_setstacksize(&tattr, 128 * 1024 * 1024);

  pthread_t th_calc;

  while(1){
    // receive UDP packets for F(a) and F(b)
    pthread_t th;

    pthread_create(&th, NULL, recv_thread, &e);
    pthread_join(th, NULL);
    swap_buffer(&e);

    // do FFT operation
    printf("kick calc\n");

    g_is_calc_finished = false;

    pthread_create(&th_calc, &tattr, calc_thread, &arg);
    pthread_detach(th_calc);

    printf("file_count = %d\n", sys.file_count);
    printf("result_count = %d\n", sys.result_count);

    if(sys.file_count * param_accumulation_time + sys.result_count >=
       param_duration){
      break;
    }
  }

  while(g_is_calc_finished == false){
    sleep(1);
  }

  env_free(&e);
  sys_free(&sys);

  return NULL;
}
