#include <mesh.h>


int main(int argc, char *argv[])
{ 
  int    iTri, iVer; int NbrEdgBoudry;
  double to, ti;
  
  if ( argc < 2 ) {
    printf(" usage : mesh file \n");
    return 0;
  }
  
  //--- read a mesh 
  to =  GetWallClock();
  int readEfr = 1;
  Mesh * Msh = msh_read(argv[1], readEfr);
  ti =  GetWallClock();
  
  if ( ! Msh ) return 0;
  
  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n",ti-to);
  //write_Efr_to_txt("Efr.txt", Msh);
   
  //--- create neigbhors Q2 version 
  to =  GetWallClock();
  msh_neighborsQ2(Msh);
  ti =  GetWallClock();
  write_TriVoi_to_file("../output/TriVoi_Q2.txt", Msh);
  printf("[Output File] neighbors written in TriVoi_Q2.txt \n");
  if (Msh->NbrEfr != 0) {   // we have boundary edge number in mesh file
    NbrEdgBoudry = Msh->NbrEfr;
  }
  else{                     // if we do not, we compute manually
    NbrEdgBoudry = compute_NbrEdgBoudry(Msh);
  }
  printf("  time q2 neigh.        %10f (s) \n",ti-to);
  printf("  Nbr boundary edges q2 %10d \n", NbrEdgBoudry);
  
  // reinialize the mesh
  free(Msh);
  Msh = msh_read(argv[1], 1);
  if ( ! Msh ) return 0;

  
  //--- create neigbhors with hash table 
  const char* keyMode = "divide";      // "min" or "sum" or "divide"
  HashTable *hsh;
  to =  GetWallClock();
  hsh = msh_neighbors(Msh, keyMode);
  ti =  GetWallClock();
  int MaxCol; double AveCol;
  collision(hsh, &MaxCol, &AveCol);
  printf("  Max collision %10d \n", MaxCol);
  printf("  Average collision %10f \n", AveCol);

  write_TriVoi_to_file("../output/TriVoi_Hash.txt", Msh);
  printf("[Output File] neighbors written in TriVoi_Hash.txt \n");
  if (Msh->NbrEfr != 0) {   // we have boundary edge number in mesh file
    NbrEdgBoudry = Msh->NbrEfr;
  }
  else{                     // if we do not, we compute manually
    NbrEdgBoudry = compute_NbrEdgBoudry(Msh);
  }
  printf("  time hash tab neigh.  %10f (s) \n",ti-to);
  printf("  Nbr boundary edges hash %10d \n", NbrEdgBoudry);
  
  //--- find connex components
  find_connex_components(Msh);

  //--- TODO: compute mesh quality
  double  *Qal1 = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));
  double  *Qal2 = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));

  msh_quality(Msh, Qal1, 1);
  msh_quality(Msh, Qal2, 2);
  
  msh_write2dfield_Triangles("../output/quality1.solb", Msh->NbrTri, Qal1);
  msh_write2dfield_Triangles("../output/quality2.solb", Msh->NbrTri, Qal2);
  printf("[Output File] quality written in quality.solb \n");
  
  //--- TODO: compute metric field
  double3d *Met = (double3d *)malloc(sizeof(double3d) * (Msh->NbrVer+1));

  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
  	 Met[iVer][0] = 1.;
  	 Met[iVer][1] = 0.;
  	 Met[iVer][2] = 1.;
  } 
  
  msh_write2dmetric("../output/metric.solb", Msh->NbrVer, Met);
  printf("[Output File] metric written in metric.solb \n");
  	
  	
  return 0;
}





