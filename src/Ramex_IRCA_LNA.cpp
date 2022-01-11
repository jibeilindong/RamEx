// Updated at Nov. 9, 2021
// Updated by Gongchao Jing
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include "utility.h" 


using namespace std;

// Parameters Def



string Infilename;
string Metafilename;
string LocalBandsAnnfile;
string Out_path = "Results_IRCA/";
string R_path;
string Error_file; 


    
int printhelp(){
    
    cout << "RamEX version: " << Version << endl;
    cout << "\tRamEX LNA: " << endl;
    cout << "Usage:" << endl;
    cout << "RamEX-IRCA-LNA [Option]" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t -i Input Raman data file" << endl;
	cout << "\t -m Input Meta-data file" << endl;
    cout << "\t -l Input the local bands annotation file" << endl;
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output path, default is \"Results_IRCA\"" << endl;
    
    cout << "\t  -h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){
    
    if (argc ==1)
	printhelp();
    R_path = Check_Env() + "/Rscript";
    LocalBandsAnnfile = Check_Env() + "/databases/" + "Local_bands_annotation.txt";
    int i = 1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            
                            case 'i':
                                 
                                 Infilename = argv[i+1];                                     
                                 break;
                                 
                            case 'm':
				 Metafilename = argv[i+1];                                     
				 break;
			    case 'l':
                                  LocalBandsAnnfile = argv[i+1];
                                  break;
                                                      
                            case 'o': Out_path = argv[i+1]; break;
                            
                 
                            case 'h': printhelp(); break;

                            default : printf("Error: Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }
    
    string command_mkdir = "mkdir " + Out_path;
    system(command_mkdir.c_str());

    Error_file = Out_path + "/error.log";
    return 0;
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    char command[BUFFER_SIZE]; 
    sprintf(command, "Rscript %s/IRCA_lna.R  -i %s -m %s -l %s -o %s", R_path.c_str(), Infilename.c_str(), Metafilename.c_str(), LocalBandsAnnfile.c_str(), Out_path.c_str());
    Run_With_Error(command, Error_file.c_str());
    cout << endl << "IRCA Local Network Analysis Finished"<< endl;

    return 0;
    }
