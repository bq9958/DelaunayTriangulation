#include <mesh.h>


int tri2edg[3][2] = { {1,2} , {2,0} , {0,1} };


Mesh * msh_init()
{
  Mesh *Msh = malloc(sizeof(Mesh));
  if ( ! Msh ) return NULL;
  
  Msh->Dim    = 0;
  Msh->NbrVer = 0;
  Msh->NbrTri = 0;
  Msh->NbrEfr = 0;
  Msh->NbrEdg = 0;

  Msh->Box[0] =  1.e30;   // xmin
  Msh->Box[1] = -1.e30;   // xmax
  Msh->Box[2] =  1.e30;   // ymin
  Msh->Box[3] = -1.e30;   // ymax
  
  //--- Data for the list of vertices
  Msh->Crd    = NULL;
  
  //--- Data for the list of triangles
  Msh->Tri    = NULL;
  Msh->TriVoi = NULL;
  Msh->TriRef = NULL;
  
  //--- Data for the list of boundary edges  
  Msh->Efr    = NULL;
  Msh->EfrVoi = NULL;
  Msh->EfrRef = NULL;

  //--- Data for the list of edges  
  Msh->Edg    = NULL;
  
  return Msh;
}  
  


Mesh * msh_read(char *file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref; 
  
  int fmsh = 0;
  
  if ( ! file ) return NULL;
  
  Mesh * Msh = msh_init();
    
  //--- set file name 
  strcpy(InpFil,file);
  if ( strstr(InpFil,".mesh") ) {
    if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".meshb");
    if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".mesh");
      if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);
  
  Msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  Msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);
  
  //--- allocate arrays   
  Msh->Crd    = calloc( (Msh->NbrVer+1), sizeof(double3d) );
  Msh->Tri    = calloc( (Msh->NbrTri+1), sizeof(int3d) );  
  Msh->TriRef = calloc( (Msh->NbrTri+1), sizeof(int1d) );  
  
  
  //--- read vertices   
  GmfGotoKwd(fmsh, GmfVertices);
  if ( Msh->Dim == 2 ) {
    if ( FilVer == GmfFloat ) {		// read 32 bits float
      for (i=1; i<=Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        Msh->Crd[i][0] = (double)bufFlt[0];
        Msh->Crd[i][1] = (double)bufFlt[1];
      }
    }
    else  {	// read 64 bits float
      for (i=1; i<=Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        Msh->Crd[i][0] = bufDbl[0];
        Msh->Crd[i][1] = bufDbl[1];
      }  
    }
  }
  else {
    fprintf(stderr,"  ## ERROR: 3D is not implemented\n");
    exit(1);
  }
  
  
  //--- read triangles   
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i=1; i<=Msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    Msh->Tri[i][0] = bufTri[0];
    Msh->Tri[i][1] = bufTri[1];
    Msh->Tri[i][2] = bufTri[2];
    Msh->TriRef[i] = bufTri[3];
  }
  
  
  //--- read boundary edges
  if ( readEfr == 1 ) {
    Msh->NbrEfr = GmfStatKwd(fmsh, GmfEdges);
    Msh->Efr    = calloc( (Msh->NbrEfr+1), sizeof(int2d) );  
    Msh->EfrRef = calloc( (Msh->NbrEfr+1), sizeof(int1d) );  

    GmfGotoKwd(fmsh, GmfEdges);
    for (i=1; i<=Msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      Msh->Efr[i][0] = bufEfr[0];
      Msh->Efr[i][1] = bufEfr[1];
      Msh->EfrRef[i] = bufEfr[2];
    }
  }
  
  
  GmfCloseMesh(fmsh);
  
  return Msh;
  
}



double * sol_read(char *file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[ GmfMaxTyp ]; 
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;
  
  int fsol = 0;
  
  if ( ! file ) return NULL;
  
  double * sol = NULL;
    
    
  //--- set file name 
  strcpy(InpFil, file);
  if ( strstr(InpFil,".sol") ) {
    if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".solb");
    if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".sol");
      if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);
  
  SolTyp = GmfSolAtVertices;		// read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);
	
	
  if ( nbrSol == 0 ) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if ( dim != mshDim ) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if ( nbrSol != mshNbrSol ) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (  NbrTyp != 1 ) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if ( TypTab[0] != GmfSca ) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }
	
	sol = (double *)calloc(nbrSol+1, sizeof(double));
	
	
  GmfGotoKwd(fsol, SolTyp);

  for (i=1; i<=nbrSol; ++i) {
		if ( FilVer == GmfFloat ) {
	    GmfGetLin(fsol, SolTyp, &bufFlt);
	    sol[i] = (double)bufFlt;
		}
		else {
	    GmfGetLin(fsol, SolTyp, &bufDbl);
	    sol[i] = bufDbl;
		}
  }
	
  if ( !GmfCloseMesh(fsol) ) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    //myexit(1);
  }
	
  return sol;	
}





int msh_boundingbox(Mesh *Msh)
{
  int1d iVer;
  
  //--- compute bounding box 
  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
    // TODO: Set Msh->Box  
  }
  
  
  return 1;
}






int msh_write(Mesh *Msh, char *file)
{
  int iVer, iTri, iEfr;
  int FilVer = 2;
  
  if ( ! Msh  ) return 0;
  if ( ! file ) return 0;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer=1; iVer<=Msh->NbrVer; iVer++) 
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0); 
  
  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++)  
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);
  
  if ( Msh->NbrEfr > 0 ) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr=1; iEfr<=Msh->NbrEfr; iEfr++)  
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }
  
  GmfCloseMesh(fmsh);
  
  return 1;
}





int msh_neighborsQ2(Mesh *Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;
  
  if ( ! Msh ) return 0;
  
  if ( Msh->TriVoi == NULL )
    Msh->TriVoi = calloc( (Msh->NbrTri+1), sizeof(int3d) );  
  
  //--- Compute the neighbors using a quadratic-complexity algorithm 
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {
      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];
			
      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices 
      for (jTri=1; jTri<=Msh->NbrTri; jTri++) {
        if ( iTri == jTri ) 
          continue;
        
        for (jEdg=0; jEdg<3; jEdg++) {
          jVer1 = Msh->Tri[jTri][ tri2edg[jEdg][0] ];
          jVer2 = Msh->Tri[jTri][ tri2edg[jEdg][1] ];
          
          // TODO: compare the 4 points 
          if ((iVer1 == iVer2 && jVer1 == jVer2) || (iVer1 == jVer2 && iVer2 == jVer1)) {
          //       set the neighbors Msh->TriVoi if both edges match
          Msh->TriVoi[iTri][iEdg] = jTri;
          Msh->TriVoi[jTri][jEdg] = iTri;
          }
        }
      }
      
    }
  }
  
  return 1;
}   



int msh_neighbors(Mesh *Msh)
{
  int iTri, iEdg, iVer1, iVer2;
  
  if ( ! Msh ) return 0;
  
  if ( Msh->TriVoi == NULL )
    Msh->TriVoi = calloc( (Msh->NbrTri+1), sizeof(int3d) );  
  
  
  //--- initialize HashTable and set the hash table 
  
  // TODO
  
  
  //--- Compute the neighbors using the hash table
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {
      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];

      // TODO:
      // compute the key : iVer1+iVer2   
      // do we have objects as that key   hash_find () */
      //  if yes ===> look among objects and potentially update TriVoi */
      //  if no  ===> add to hash table   hash_add()   */
    }
  }
  
  return 1;
}   

int msh_quality(Mesh *Msh, double *Qal, int mode)
{
  double alpha1 = 4.0/sqrt(3.0); double alpha2 = 6.0/sqrt(3.0);
  double l1, l2, l3, area;

  for (int iTri=1; iTri<=Msh->NbrTri; iTri++) {
    int iVer0 = Msh->Tri[iTri][0]; int iVer1 = Msh->Tri[iTri][1]; int iVer2 = Msh->Tri[iTri][2];
    double x0 = Msh->Crd[iVer0][0]; double y0 = Msh->Crd[iVer0][1];
    double x1 = Msh->Crd[iVer1][0]; double y1 = Msh->Crd[iVer1][1];
    double x2 = Msh->Crd[iVer2][0]; double y2 = Msh->Crd[iVer2][1];
    
    l1 = distance(x0, y0, x1, y1); l2 = distance(x1, y1, x2, y2); l3 = distance(x2, y2, x0, y0);
    area = triArea(x0, y0, x1, y1, x2, y2); 

    if (mode == 1){
      Qal[iTri] = alpha1 * (pow(l1,2) + pow(l2,2) + pow(l3,2)) / area;
    }
    else if (mode == 2){
      double rho = 2.0 * area / (l1 + l2 + l3);
      Qal[iTri] = alpha2 * fmax(fmax(l1, l2), l3) / rho;
    }
  } 

 return 0;
}


HashTable * hash_init(int SizHead, int NbrMaxObj)
{
	HashTable *hsh = NULL;
	
	// to be implemented

	// allocate hash table
	
	// initialize hash table
	
	// allocate Head, LstObj
	
	
  return hsh;
}


int hash_find(HashTable *hsh, int iVer1, int iVer2)
{
  
	// to be implemented
	
	// return the id found (in LstObj ), if 0 the object is not in the list
	
	return 0;
}


int hash_add(HashTable *hsh, int iVer1, int iVer2, int iTri)
{

  // to be implemented
	
  // ===> add this entry in the hash tab 
	
	return 0;
}




int hash_suppr(HashTable *hsh, int iVer1, int iVer2, int iTri)
{

  // to be implemented
	
  // ===> suppress this entry in the hash tab 
	
	return 0;
}


int msh_write2dfield_Vertices(char *file, int nfield, double *field) 
{
  int iVer;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSca;
  
  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);
  
  for (iVer=1; iVer<=nfield; iVer++) 
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}



int msh_write2dfield_Triangles(char *file, int nfield, double *field) 
{
  int iTri;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSca;
  
  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);
  
  for (iTri=1; iTri<=nfield; iTri++) 
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}



int msh_write2dmetric(char *file, int nmetric, double3d *metric) 
{  
  int iVer;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSymMat;
  
  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);
  
  for (iVer=1; iVer<=nmetric; iVer++) 
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double triArea(double x0, double y0, double x1, double y1, double x2, double y2)
{
    return 0.5 * ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
}


