#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>


using namespace std; 

int main ()
{
  int cptfile = 1;
  int result;
  vector<string> tache = {"30","60"};
  vector<string> scale = {"1","10"};
  vector<string> k = {"1_25","2_5","5_0"};
  vector<string> df = {"0_1","0_15","0_2","0_25"};
  vector<string> mf = {"0_1","0_2","0_4","0_8","1_6","3_2"};

  for(const auto& t : tache){
    for(const auto& sc : scale){
        for(const auto& kk : k){
            for(const auto& d : df){
                for(const auto& m : mf){
                    string oldname = "inst"+t+"-"+kk+"-"+d+"-"+m+"-"+sc+".dat";
                    string newname = "inst_"+to_string(cptfile);
                    result = rename(oldname.c_str(),newname.c_str());
                    if(result != 0) cout << "erreur "<<oldname<<endl;
                    ++cptfile;
                }
            }
        }
    }
  }
  return 0;
}