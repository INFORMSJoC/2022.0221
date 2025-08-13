// basic file operations
#include <iostream>
#include <fstream>
#include <limits>       // std::numeric_limits

using namespace std;

int main (int argc, char *argv[]) {
  ifstream myfile;
  myfile.open (argv[1]);
  int nbOfTasks;
 
  myfile >> nbOfTasks;
  int *dem = new int[nbOfTasks];
  int *dur = new int[nbOfTasks];
  int *r = new int[nbOfTasks];
  int *d = new int[nbOfTasks];
  int mindem=std::numeric_limits<int>::max();
  int maxdem=0;
  int sumdem=0;
  int maxd=0;
  for (int i=0;i<nbOfTasks;i++) {
	myfile >> dem[i] >> dur[i] >> r[i] >> d[i] ;
        if (mindem > dem[i])
           mindem=dem[i];
        if (maxdem < dem[i])
           maxdem=dem[i];
        sumdem+=dem[i];
        if (maxd < d[i])
           maxd=d[i];
  }
  myfile.close();
  int maxmaxdemperperiod=0;
  int minmaxdemperperiod=sumdem;
  int totmaxdemperperiod = 0;
  for (int t=0;t<maxd;t++) {
     int sumdem_t=0;
     for (int i=0;i<nbOfTasks;i++)
        if (t>=r[i] && t<d[i])
           sumdem_t+=dem[i];
     totmaxdemperperiod = totmaxdemperperiod + sumdem_t;
     if (maxmaxdemperperiod < sumdem_t)
        maxmaxdemperperiod=sumdem_t;
    if (minmaxdemperperiod > sumdem_t)
        minmaxdemperperiod=sumdem_t;
  }
  
    cout << argv[1]  << "\t horizon=" << maxd << "\n per task : \t mindem=" << mindem << "\t maxdem=" << maxdem << "\t totdem=" << sumdem << "\t avgdem=" << (float) sumdem / (float) nbOfTasks << "\n per period: \t minmaxdemperperiod=" << minmaxdemperperiod << "\t maxmaxdemperperiod=" << maxmaxdemperperiod << "\t totmaxdemperperiod=" <<  totmaxdemperperiod << "\t avgmaxdemperperiod=" << (float) totmaxdemperperiod / (float) maxd << endl;
  return 0;
}
