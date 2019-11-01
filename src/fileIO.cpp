#include "fileIO.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

int fileIO::CountNumbersOfTextLines(const std::string &filePath )
{
  long i = 0;
  std::ifstream ifs( filePath );
  if( ifs ){
    std::string line;

    while( true ){
      std::getline( ifs, line );
      i++;
      if( ifs.eof() )
      break;
    }
  }

  return i-1;
}

void fileIO::read_original_file_double(DOUBLEARRAY2 &x,const int &nd,const int &number,const std::string &file)
{
  std::string str,tmp;
  std::ifstream node_file(file);
  if(!node_file){
    std::cout << "Error:Input "<< file << " not found" << std::endl;
    exit(1); 
  }
  for(int i=0;i<nd;i++){

    std::getline(node_file,str);
    std::istringstream stream(str);

    for(int j=0;j<number;j++){
      std::getline(stream,tmp,' ');
      x[i][j] = std::stod(tmp);
    }
  }
}

void fileIO::read_original_file_int(INTARRAY2 &x,const int &nd,const int &number,const std::string &file)
{
  std::string str,tmp;
  std::ifstream node_file(file);
  if(!node_file){
    std::cout << "Error:Input "<< file << " not found" << std::endl;
    exit(1); 
  }
  for(int i=0;i<nd;i++){

    std::getline(node_file,str);
    std::istringstream stream(str);

    for(int j=0;j<number;j++){
      std::getline(stream,tmp,' ');
      x[i][j] = std::stoi(tmp);
    }
  }
}

