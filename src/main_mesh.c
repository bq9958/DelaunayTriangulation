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
  double  *Qal1 = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));
  double  *Qal2 = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));

  msh_quality(Msh, Qal1, 1);
  msh_quality(Msh, Qal2, 2);
  
  msh_write2dfield_Triangles("quality1.solb", Msh->NbrTri, Qal1);
  msh_write2dfield_Triangles("quality2.solb", Msh->NbrTri, Qal2);
  printf("  quality written in quality.solb \n");
  
  //--- TODO: compute metric field
  double3d *Met = (double3d *)malloc(sizeof(double3d) * (Msh->NbrVer+1));

  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
  	 Met[iVer][0] = 1.;
  	 Met[iVer][1] = 0.;
  	 Met[iVer][2] = 1.;
  } 
  
  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);
  printf("  metric written in metric.solb \n");
  	
  	
  return 0;
}





