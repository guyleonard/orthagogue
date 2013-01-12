#include <sys/types.h> // some systems require it
#include <sys/stat.h>
#include <sys/termios.h> // for winsize
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/times.h>
#define MAXLINE 4096

bool next_proc = true; // Sets to false if a process abrupst abnormally
int n_cpus = 1;
float cpu_1_t_real = 1.0;
uint total_file_cnt = 3; // The number to iterate of if the param is '-1'
static void err_sys(char *err) {
  printf("!!\n!!\tAborts due to %s. Continues on the enxt run\n", err);
    next_proc =false;
  //  exit(2);
}

static double t_real = 0, t_user = 0, t_sys = 0, t_c_sys = 0, t_c_user=0;
int iterations = 0;
static void pr_exit(int stat) {//printf("ok\tProsess return %d as an exit status\n", stat);
  if(stat == 2) 
    next_proc =false;
  //  exit(2); 
}
long clktck = 0;
 void pr_times(clock_t real, struct tms *tmsstart, struct tms *tmsend, const bool print_data) {
  if(!print_data) {
    iterations++;

    if(clktck == 0) // henter antall tikk foerte gang
      if((clktck =  sysconf(_SC_CLK_TCK)) < 0)
	err_sys("sysconf error");
    t_real += real/(double)clktck;
    t_user += (tmsend->tms_utime - tmsstart->tms_utime) / (double) clktck;
    t_sys += (tmsend->tms_stime - tmsstart->tms_stime) / (double) clktck;
    t_c_sys += (tmsend->tms_cstime - tmsstart->tms_cstime) / (double) clktck;
    t_c_user += (tmsend->tms_cutime - tmsstart->tms_cutime) / (double) clktck;

  } else {
    if(n_cpus == 1) cpu_1_t_real = t_real;
    if(true) {
      const float sec =  t_real/iterations;
      printf("-\t%.3f%c is the parallelism using %d threads. (Total time is %.2f seconds = %.2f minutes.)\n", (float)100*(cpu_1_t_real/(t_real*n_cpus)), '%', n_cpus, sec, sec/60);
    } else {
      printf("\n# Run %d times the program, giving the following data as an average:\n", iterations);
      printf(" real: %7.4f\n", t_real/iterations);
      printf(" user: %7.4f\n", 	   t_user/iterations);
      printf(" sys: %7.4f\n", t_sys/(double)iterations);
      printf(" child user:  %7.4f\n", t_c_user/(double)iterations);
      printf(" child sys:  %7.4f\n", t_c_sys/(double)iterations);
    }
  }
}

//tmsstart = tms(), tmsend = tms();
static void do_cmd(char *cmd) {
  struct tms tmsstart, tmsend;
  clock_t start, end;
  int status;
  //  printf("\nCommano: %s\n", cmd);
  if((start = times(&tmsstart)) == -1) // starting values
    err_sys("times error");
  if((status = system(cmd)) < 0)
    err_sys("system() error");
  if((end = times(&tmsend)) == -1)
    err_sys("times error\n");
  pr_times(end-start, &tmsstart, &tmsend, false);
  pr_exit(status);
}

char *get_str(int id) {
  // must uupdate 'total_file_cnt = 3'; // The number to iterate of if the param is '-1'
  //  if(id == -1) return ".././orthAgogue -i  /work/mironov/outputs/mpiblast/AtCeDmHsMm/AtCeDmHsMm.blast -p 0 -t 1 -s '_'";
  if(id == 0) return "./orthAgogue -i  all.blast -p 0 -t 1 -s '_'";
  else if(id == 1) return "./orthAgogue -i  /work/mironov/outputs/mpiblast/AtCeDmHsMm/AtCeDmHsMm.blast -p 0 -t 1 -s '_'";
  else if(id == 2) return "./orthAgogue -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' ";
  else {
    fprintf(stderr, "!!\tid %d given is not implemented. Aborts\n", id);
    exit(2);
  }
}
char *get_exec_term(int id, int n_cpus) {
  char *string = new char[400]; memset(string, '\0', 400);
  sprintf(string, "%s -c %d", get_str(id), n_cpus);
  return string;
}
int main(int argc, char *argv[]) {
  if(argc < 2) printf("!!\tName of program not set: aborts execution.\n");
  else {
    const int id = atoi(argv[1]);
    if(id != -1) {
      const int max_cpus = 10;
      printf("\nThe result of the data set, using at max %d cpus for '%s' is:\n", max_cpus, get_str(id));
      for(n_cpus = 1; n_cpus < max_cpus; n_cpus++) {    
	char *string = get_exec_term(id, n_cpus);    
	clktck = 0;
	t_real = 0, t_user = 0, t_sys = 0, t_c_sys = 0, t_c_user=0;
	iterations = 0;
	setbuf(stdout, NULL);
	for(int i = 0; i < 3; i++) { // Antall ganger programmet skal kjore
	  if(next_proc) do_cmd(string); // Bruker commando-linjen
	}
	delete [] string;
	pr_times(0, NULL, NULL, true);
      }
    } else {
      const int max_cpus = 10;
      printf("total_file_cnt = %d\n", total_file_cnt);
      for(int id = 0; id < (int)total_file_cnt; id++) {
	printf("The result of the data set, using %d cpus for '%s' is:\n", max_cpus, get_str(id));
	for(n_cpus = 1; n_cpus < max_cpus; n_cpus++) {    
	  char *string = get_exec_term(id, n_cpus);    
	  clktck = 0;
	  t_real = 0, t_user = 0, t_sys = 0, t_c_sys = 0, t_c_user=0;
	  iterations = 0;
	  setbuf(stdout, NULL);
	  for(int i = 0; i < 3; i++) { // Antall ganger programmet skal kjore
	    if(next_proc) do_cmd(string); // Bruker commando-linjen
	  }
	  delete [] string;
	  pr_times(0, NULL, NULL, true);
	}
      }
    }
    exit(0);
  }
  return 0;
}
