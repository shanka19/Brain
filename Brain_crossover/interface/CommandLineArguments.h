#include "string"
#include "map"

std::map<std::string, float> commandLineArguments(int argc, char *argv[])
{
  std::map<std::string, float> cmdMap;
  for (unsigned int i=1; i<argc; i+=2) 
  {
    cmdMap[std::string(argv[i])]=atof(argv[i+1]);
  }
  return cmdMap;
}
