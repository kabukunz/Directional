#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/line_cylinders.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/branched_isolines.h>

#include <directional/mesh_function_isolines.h>

#include <directional/setup_mesh_function_isolines.h>
#include "polygonal_write_OFF.h"

#define NUM_N 1

int N[NUM_N];
int currN = 0;
Eigen::MatrixXi FMeshWhole, FMeshCut[NUM_N], FField[NUM_N], FSings[NUM_N], FSeams[NUM_N], FIso[NUM_N];
Eigen::MatrixXd VMeshWhole, VMeshCut[NUM_N], VField[NUM_N], VSings[NUM_N], VSeams[NUM_N], VIso[NUM_N];
Eigen::MatrixXd CField[NUM_N], CSeams[NUM_N], CSings[NUM_N], CIso[NUM_N];
Eigen::MatrixXd rawField[NUM_N], combedField[NUM_N];
// Eigen::MatrixXd rawField[NUM_N], combedField[NUM_N], barycenters;
Eigen::VectorXd effort[NUM_N], combedEffort[NUM_N];
Eigen::VectorXi matching[NUM_N], combedMatching[NUM_N];
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi DPolyMesh[NUM_N];
Eigen::MatrixXi FPolyMesh[NUM_N];
Eigen::MatrixXd VPolyMesh[NUM_N];
Eigen::VectorXi singIndices[NUM_N], singVertices[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
igl::opengl::glfw::Viewer viewer;



typedef enum {FIELD, INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_raw_field_mesh()
{
  for (int i=1;i<=3;i++)  //hide all other meshes
    viewer.data_list[i].show_faces=(viewingMode==FIELD);
  
  viewer.data_list[4].show_faces=(viewingMode==INTEGRATION);
  
  if (viewingMode==FIELD){
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField[currN], FField[currN]);
    viewer.data_list[1].set_colors(CField[currN]);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings[currN], FSings[currN]);
    viewer.data_list[2].set_colors(CSings[currN]);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    viewer.data_list[3].clear();
    viewer.data_list[3].set_mesh(VSeams[currN], FSeams[currN]);
    viewer.data_list[3].set_colors(CSeams[currN]);
    viewer.data_list[3].show_faces = true;
    viewer.data_list[3].show_lines = false;
  } else {
    viewer.data_list[4].clear();
    viewer.data_list[4].set_mesh(VIso[currN], FIso[currN]);
    viewer.data_list[4].set_colors(CIso[currN]);
    viewer.data_list[4].show_faces = true;
    viewer.data_list[4].show_lines = false;
  }
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = INTEGRATION; break;
    case '3': currN=(currN+1)%NUM_N; break;
  }
  update_raw_field_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show isoline mesh" << std::endl <<
  "  3  change between different N" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/vase.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-4.rawfield", N[0], rawField[0]);
//   directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-7.rawfield", N[1], rawField[1]);
//   directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-11.rawfield", N[2], rawField[2]);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
//   igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  bool verbose=true;
  
  //combing and cutting
  for (int i=0;i<NUM_N;i++){
    
    
    // COMBING


    // principal matching
    directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField[i], matching[i], effort[i]);

    // singularities
    directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort[i],matching[i], N[i],singVertices[i], singIndices[i]);
    

    // INTEGRATION


    // integration
    std::cout<<"Setting up Integration #"<<i<<std::endl;
    
    directional::IntegrationData intData(N[i]);
    directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField[i], matching[i], singVertices[i], intData, VMeshCut[i], FMeshCut[i], combedField[i], combedMatching[i]);
    
    intData.verbose=false;
    intData.integralSeamless=true;
    intData.roundSeams=false;
  
    std::cout<<"Solving integration for N="<<N[i]<<std::endl;
    
    directional::integrate(VMeshWhole, FMeshWhole, FE, combedField[i],  intData, VMeshCut[i], FMeshCut[i], NFunction[i],NCornerFunction[i]);
    
    std::cout<<"Done!"<<std::endl;


    // MESHING


    // The meshing unit is independent from the integration unit, and can be potentially used with external functions; one should fill the MeshFunctionIsolinesData structure with the input, and call mesh_function_isolines().    
    
    //setting up mesh data from itnegration data
    directional::MeshFunctionIsolinesData mfiData;
    directional::setup_mesh_function_isolines(VMeshCut[i], FMeshCut[i], intData, mfiData);
    
    //meshing
    directional::mesh_function_isolines(VMeshWhole, FMeshWhole,EV, EF, FE, mfiData,  verbose, VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);
    

    
    // SAVING


    //saving
    hedra::polygonal_write_OFF(TUTORIAL_SHARED_PATH "/vase-"+std::to_string(N[i])+"-generated.off", VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);
    
    
    // SHOWING

    //raw field mesh
    directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField[i], directional::indexed_glyph_colors(combedField[i]), VField[i], FField[i], CField[i],1.0);
    
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[1].clear();
      viewer.data_list[1].set_mesh(VField[i], FField[i]);
      viewer.data_list[1].set_colors(CField[i]);
      viewer.data_list[1].show_faces = true;
      viewer.data_list[1].show_lines = false;
    }
    
    //singularity mesh
    directional::singularity_spheres(VMeshWhole, FMeshWhole, N[i], singVertices[i], singIndices[i], VSings[i], FSings[i], CSings[i],2.5);
    
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[2].clear();
      viewer.data_list[2].set_mesh(VSings[i], FSings[i]);
      viewer.data_list[2].set_colors(CSings[i]);
      viewer.data_list[2].show_faces = true;
      viewer.data_list[2].show_lines = false;
    }
    
    //seams mesh
    Eigen::VectorXi isSeam=Eigen::VectorXi::Zero(EV.rows());
    for (int i=0;i<FE.rows();i++)
      for (int j=0;j<3;j++)
        if (intData.face2cut(i,j))
          isSeam(FE(i,j))=1;
    directional::seam_lines(VMeshWhole, FMeshWhole, EV, combedMatching[i], VSeams[i], FSeams[i], CSeams[i],2.5);
    
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[3].clear();
      viewer.data_list[3].set_mesh(VSeams[i], FSeams[i]);
      viewer.data_list[3].set_colors(CSeams[i]);
      viewer.data_list[3].show_faces = true;
      viewer.data_list[3].show_lines = false;
    }
    
    directional::branched_isolines(VMeshCut[i], FMeshCut[i],NFunction[i], VIso[i], FIso[i], CIso[i]);
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[4].clear();
      viewer.data_list[4].set_mesh(VIso[i], FIso[i]);
      viewer.data_list[4].set_colors(CIso[i]);
      viewer.data_list[4].set_face_based(true);
      viewer.data_list[4].show_faces = false;
      viewer.data_list[4].show_lines = false;
    }
  }
  
  
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
  viewer.data_list[0].set_colors(directional::default_mesh_color());
  viewer.data_list[0].set_face_based(false);
  viewer.data_list[0].show_lines=false;
  update_raw_field_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


