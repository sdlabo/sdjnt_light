#include <stdio.h>
#include <signal.h>
#include <string>
#include "sdlab_signal.h"
#include "cmdline.h"

using namespace std;

char param_basename[1024];
int param_duration;
int param_accumulation_time;

void sigint_handler(int sig)
{
  exit(0);
}

int main(int argc, char **argv)
{
  cmdline::parser p;
  p.add<string>("basename", 'b', "base name of filename. "
                "Output filename will be [basename]_[count].",
                false,
                "sdjnt_light_default");
  p.add<int>("duration", 'd', "measurement duration [sec].",
             false, 10);
  p.add<int>("accumulation", 'a',
             "the number of accumulated results per file.",
             false, 5);
  p.add("help", 'h', "print help");

  p.parse(argc, argv);
  if(p.parse(argc, argv) == false || p.exist("help")){
    printf("%s\n", p.error_full().c_str());
    printf("%s", p.usage().c_str());
    return 0;
  }


  param_duration = p.get<int>("duration");
  param_accumulation_time = p.get<int>("accumulation");
  strcpy(param_basename, p.get<string>("basename").c_str());

  printf("param_dulation = %d\n", param_duration);
  printf("param_accumulation_time = %d\n", param_accumulation_time);
  printf("param_basename = %s\n", param_basename);

  signal(SIGINT, sigint_handler);
  sdlab_signal_main();

  return 0;
}
