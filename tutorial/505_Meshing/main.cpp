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

int N;
int currN = 0;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams, FIso;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams, VIso;
Eigen::MatrixXd CField, CSeams, CSings, CIso;
Eigen::MatrixXd rawField, combedField;
// Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi DPolyMesh;
Eigen::MatrixXi FPolyMesh;
Eigen::MatrixXd VPolyMesh;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd NFunction, NCornerFunction;
igl::opengl::glfw::Viewer viewer;



typedef enum {FIELD, INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_raw_field_mesh()
{
  for (int i=1;i<=3;i++)  //hide all other meshes
    viewer.data_list.show_faces=(viewingMode==FIELD);
  
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
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-4.rawfield", N, rawField);
//   directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-7.rawfield", N[1], rawField[1]);
//   directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-11.rawfield", N[2], rawField[2]);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
//   igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  bool verbose=true;
  
  //combing and cutting
  for (int i=0;i<NUM_N;i++){
    
    
    // COMBING


    // principal matching
    directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort);

    // singularities
    directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,singVertices, singIndices);
    

    // INTEGRATION


    // integration
    std::cout<<"Setting up Integration #"<<i<<std::endl;
    
    directional::IntegrationData intData;
    directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
    
    intData.verbose=false;
    intData.integralSeamless=true;
    intData.roundSeams=false;
  
    std::cout<<"Solving integration for N="<<N<<std::endl;
    
    directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut, NFunction,NCornerFunction);
    
    std::cout<<"Done!"<<std::endl;


    // MESHING


    // The meshing unit is independent from the integration unit, and can be potentially used with external functions; one should fill the MeshFunctionIsolinesData structure with the input, and call mesh_function_isolines().    
    
    //setting up mesh data from itnegration data
    directional::MeshFunctionIsolinesData mfiData;
    directional::setup_mesh_function_isolines(VMeshCut, FMeshCut, intData, mfiData);
    
    //meshing
    directional::mesh_function_isolines(VMeshWhole, FMeshWhole,EV, EF, FE, mfiData,  verbose, VPolyMesh, DPolyMesh, FPolyMesh);
    

    
    // SAVING


    //saving
    hedra::polygonal_write_OFF(TUTORIAL_SHARED_PATH "/vase-4-generated.off", VPolyMesh, DPolyMesh, FPolyMesh);
    
    
    // SHOWING

    //raw field mesh
    directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, directional::indexed_glyph_colors(combedField), VField, FField, CField,1.0);
    
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[1].clear();
      viewer.data_list[1].set_mesh(VField, FField);
      viewer.data_list[1].set_colors(CField);
      viewer.data_list[1].show_faces = true;
      viewer.data_list[1].show_lines = false;
    }
    
    //singularity mesh
    directional::singularity_spheres(VMeshWhole, FMeshWhole, N, singVertices, singIndices, VSings, FSings, CSings,2.5);
    
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[2].clear();
      viewer.data_list[2].set_mesh(VSings, FSings);
      viewer.data_list[2].set_colors(CSings);
      viewer.data_list[2].show_faces = true;
      viewer.data_list[2].show_lines = false;
    }
    
    //seams mesh
    Eigen::VectorXi isSeam=Eigen::VectorXi::Zero(EV.rows());
    for (int i=0;i<FE.rows();i++)
      for (int j=0;j<3;j++)
        if (intData.face2cut(i,j))
          isSeam(FE(i,j))=1;
    directional::seam_lines(VMeshWhole, FMeshWhole, EV, combedMatching, VSeams, FSeams, CSeams,2.5);
    
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[3].clear();
      viewer.data_list[3].set_mesh(VSeams, FSeams);
      viewer.data_list[3].set_colors(CSeams);
      viewer.data_list[3].show_faces = true;
      viewer.data_list[3].show_lines = false;
    }
    
    directional::branched_isolines(VMeshCut, FMeshCut,NFunction, VIso, FIso, CIso);
    if (i==0){
      viewer.append_mesh();
      viewer.data_list[4].clear();
      viewer.data_list[4].set_mesh(VIso, FIso);
      viewer.data_list[4].set_colors(CIso);
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


