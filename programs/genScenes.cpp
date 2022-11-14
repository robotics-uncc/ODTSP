#include <stdio.h>
#include <string>
#include <iostream>
#include <MathTools.h>

int main(int argc, char** argv) {
  printf("--- Scene generator for ODTSP --- \n");
  printf("How to use: \n$ genScenes : follow on-screen instructions. \n");
  printf("$ genScenes prefix numTargets numScenes : use default settings \n\n");
  
  
  // enter default values
  double boxDim = 15.0;
  int numTargs = 10;
  int precision = 1000;
  double orbitAngle = 180.0; 
  double aspectAngleMaxDeg = 20.0;
  double orbitRadius = 1.0; 
  double orbitRadiusMax = 3.0;
  double minTurnRadius = 0.5;
  double windDirDeg = 90.0;
  double windMagnitude = 0.2;
  double vehicleSpeed = 1.0;

  
  // ask user 
  int defaultFlag = 0;
  char prefix[100];
  std::string prefixString;
  int numScenes;
  if ( argc == 1 ){
    std::cout << "Use default parameters (0-no,1-yes): "; std::cin  >> defaultFlag;
    if ( defaultFlag == 0 ){
      std::cout << "Enter the dimensions of the box : "; std::cin  >> boxDim; 
      std::cout << "Enter the orbit angle (deg) : "; std::cin  >> orbitAngle; 
      std::cout << "Enter the min. orbit radius (m) : "; std::cin  >> orbitRadius; 
      std::cout << "Enter the max. orbit radius (m) : "; std::cin  >> orbitRadiusMax; 
      std::cout << "Enter the max. pointing error (deg) : "; std::cin  >> aspectAngleMaxDeg;
      std::cout << "Enter the turn radius (m) : "; std::cin  >> minTurnRadius; 
      std::cout << "Enter the vehicle speed (m/s) : "; std::cin  >> vehicleSpeed; 
      std::cout << "Enter the wind speed (m/s) : "; std::cin  >> windMagnitude; 
      std::cout << "Enter the wind direction (deg) : "; std::cin  >> windDirDeg; 
      std::cout << "Enter the precision : "; std::cin  >> precision; 
    }    
    std::cout << "Enter the prefix to be used for scenes : "; std::cin  >> prefix;  
    std::cout << "Enter the no. of targets : "; std::cin  >> numTargs;
    std::cout << "Enter number of scenes to generate : "; std::cin  >> numScenes; 
  }
  else if ( argc == 4 ){
    prefixString = argv[1];
    strcpy(prefix, prefixString.c_str());
    numTargs = atoi(argv[2]);
    numScenes = atoi(argv[3]);    
  }
  else {
    printf("Error with input arguments -- Aborting.\n");
    return 0;
  }
  double xmin = 0; double xmax = boxDim;
  double ymin = 0; double ymax = boxDim;
  for (int i = 0 ; i < numScenes; i++ ){
    char buffer[100];
    sprintf(buffer,"%s_%04d.scn",prefix,i);
    std::string filename = buffer;
    FILE * fid = fopen(filename.c_str(),"w");
    fprintf(fid,"FILENAME\n%s\n",filename.c_str());
    fprintf(fid,"ORBIT_ANGLE_DEG\n%3.3f\n",orbitAngle);
    fprintf(fid,"ORBIT_RADIUS_MIN\n%3.3f\n",orbitRadius);
    fprintf(fid,"ORBIT_RADIUS_MAX\n%3.3f\n",orbitRadiusMax);
    fprintf(fid,"POINTING_ERROR_MAX_DEG\n%3.3f\n",aspectAngleMaxDeg);
    fprintf(fid,"MIN_TURN_RADIUS\n%3.3f\n",minTurnRadius);
    fprintf(fid,"VEHICLE_SPEED\n%3.3f\n",vehicleSpeed);
    fprintf(fid,"WIND_SPEED\n%3.3f\n",windMagnitude);
    fprintf(fid,"WIND_DIR_DEG\n%3.3f\n",windDirDeg);
    fprintf(fid,"PRECISION\n%d\n",precision);
    fprintf(fid,"BOX_DIM\n%3.3f\n",boxDim);  
    fprintf(fid,"NUM_TARGS\n%d\n",numTargs);
    fprintf(fid,"TARGS\n");
    std::vector<double> xTargs, yTargs;
    if ( numTargs > 0 || boxDim <= 0 ){
      xTargs = MathTools::randomNumberUniformDistribution( xmin, xmax, numTargs); 
      yTargs = MathTools::randomNumberUniformDistribution( ymin, ymax, numTargs);
    }
    else{
      printf("Invalid entry.\n");
      return 0;
    }
    for (int j = 0; j < numTargs; j++){
      fprintf(fid,"%3.4f %3.4f\n",xTargs[j],yTargs[j]);
    }
    fprintf(fid,"EOF\n");
    fclose(fid);
  }


}
