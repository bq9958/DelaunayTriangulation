#include <mesh.h>


int main(int argc, char *argv[])
{ 
  int    iTri, iVer;
  double to, ti;
  
  if ( argc < 2 ) {
    printf(" usage : mesh file \n");
    return 0;
  }
  
  //--- read a mesh 
  to =  GetWallClock();
  Mesh * Msh = msh_read(argv[1], 0);
  ti =  GetWallClock();
  
  if ( ! Msh ) return 0;
  
  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n",ti-to);
  
   
  //--- create neigbhors Q2 version 
  to =  GetWallClock();
  msh_neighborsQ2(Msh);
  ti =  GetWallClock();
  printf("  time q2 neigh.        %lg (s) \n",ti-to);
  
  
  //--- create neigbhors with hash table 
  to =  GetWallClock();
  msh_neighbors(Msh);
  ti =  GetWallClock();
  printf("  time hash tab neigh.  %lg (s) \n",ti-to);
  
  
  //--- TODO: compute mesh quality
  double  *Qal = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));
  
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
  	 Qal[iTri] = (double)iTri / 10.;
  } 
  
  msh_write2dfield_Triangles("quality.solb", Msh->NbrTri, Qal);

  
  //--- TODO: compute metric field
  double3d *Met = (double3d *)malloc(sizeof(double3d) * (Msh->NbrVer+1));

  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
  	 Met[iVer][0] = 1.;
  	 Met[iVer][1] = 0.;
  	 Met[iVer][2] = 1.;
  } 
  
  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);
  	
  	
  return 0;
}





