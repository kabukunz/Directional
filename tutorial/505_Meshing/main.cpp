#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <igl/Timer.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/line_cylinders.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/branched_isolines.h>

#include <directional/mesh_function_isolines.h>

#include <directional/setup_mesh_function_isolines.h>
#include <directional/directional_viewer.h>
#include "polygonal_write_OFF.h"

#define NUM_N 1

int N;
int currN = 0;
Eigen::MatrixXi FMeshWhole, FMeshCut[NUM_N];
Eigen::MatrixXd VMeshWhole, VMeshCut[NUM_N];
Eigen::MatrixXd rawField[NUM_N], combedField[NUM_N];
Eigen::VectorXd effort[NUM_N], combedEffort[NUM_N];
Eigen::VectorXi matching[NUM_N], combedMatching[NUM_N];
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi DPolyMesh[NUM_N];
Eigen::MatrixXi FPolyMesh[NUM_N];
Eigen::MatrixXd VPolyMesh[NUM_N];
Eigen::VectorXi singIndices[NUM_N], singVertices[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
directional::DirectionalViewer viewer;

typedef enum {FIELD, INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_viewer()
{
  for (int i=0;i<NUM_N;i++){
    viewer.toggle_field(false,i);
    viewer.toggle_singularities(false,i);
    viewer.toggle_seams(false,i);
    viewer.toggle_isolines(false,i);
  }
  if (viewingMode==FIELD){
    viewer.toggle_field(true,currN);
    viewer.toggle_singularities(true,currN);
    viewer.toggle_seams(true,currN);
    viewer.toggle_isolines(false,currN);
  } else {
    viewer.toggle_field(false,currN);
    viewer.toggle_singularities(false,currN);
    viewer.toggle_seams(false,currN);
    viewer.toggle_isolines(true,currN);
  }
}

    viewer.data_list[4].show_faces = (viewingMode == INTEGRATION);

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
  update_viewer();
  return true;
}

        viewer.data_list[2].clear();
        viewer.data_list[2].set_mesh(VSings, FSings);
        viewer.data_list[2].set_colors(CSings);
        viewer.data_list[2].show_faces = true;
        viewer.data_list[2].show_lines = false;

int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show isoline mesh" << std::endl <<
  "  3  change between different N" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/vase.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-4.rawfield", N[0], rawField[0]);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-7.rawfield", N[1], rawField[1]);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-11.rawfield", N[2], rawField[2]);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  
  bool verbose=true;
  
  //combing and cutting
  for (int i=0;i<NUM_N;i++){
    directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField[i], matching[i], effort[i],singVertices[i], singIndices[i]);

    directional::IntegrationData intData(N[i]);
    std::cout<<"Setting up Integration #"<<i<<std::endl;
    directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField[i], matching[i], singVertices[i], intData, VMeshCut[i], FMeshCut[i], combedField[i], combedMatching[i]);
    
    intData.verbose=false;
    intData.integralSeamless=true;
    intData.roundSeams=false;
  
    std::cout<<"Solving integration for N="<<N[i]<<std::endl;
    directional::integrate(VMeshWhole, FMeshWhole, FE, combedField[i],  intData, VMeshCut[i], FMeshCut[i], NFunction[i],NCornerFunction[i]);
    
    std::cout<<"Done!"<<std::endl;
    
    //setting up mesh data from itnegration data
    directional::MeshFunctionIsolinesData mfiData;
    directional::setup_mesh_function_isolines(VMeshCut[i], FMeshCut[i], intData, mfiData);
    
    //meshing and saving
    directional::mesh_function_isolines(VMeshWhole, FMeshWhole,EV, EF, FE, mfiData,  verbose, VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);
    hedra::polygonal_write_OFF(TUTORIAL_SHARED_PATH "/vase-"+std::to_string(N[i])+"-generated.off", VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);
    
    
    viewer.set_mesh(VMeshWhole, FMeshWhole,i);
    viewer.set_field(combedField[i], directional::DirectionalViewer::indexed_glyph_colors(combedField[i]), i);
    viewer.set_singularities(singVertices[i], singIndices[i]);
    viewer.set_seams(combedMatching[i], i);
    viewer.set_isolines(VMeshCut[i], FMeshCut[i],NFunction[i],i);
  }
  
  update_viewer();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

int main(int argc, char *argv[])
{
    // check args
    if (argc < 2)
    {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " <offfile>" << std::endl;
        return 1;
    }

    std::string cmdlineFile = argv[1];

    std::string dataDir = std::string("C:\\Users\\kabukunz\\Developer\\Content\\Data\\");
    infilename = dataDir + cmdlineFile + ".off";
    rawfieldname = dataDir + cmdlineFile + ".rawfield";
    outfilename = dataDir + cmdlineFile + "-generated.off";

    igl::readOFF(infilename, VMeshWhole, FMeshWhole);
    directional::read_raw_field(rawfieldname, N, rawField);
    igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

    std::cout << "  1  Loaded field" << std::endl
              << "  2  Show isoline mesh" << std::endl
              << "  3  change between different N" << std::endl;

    bool verbose = true;

    //combing and cutting
    for (int i = 0; i < NUM_N; i++)
    {

        // TODO: rawfield generation

        // TIMER
        timer.start();

        // COMBING

        // principal matching
        directional::principal_matching(VMeshWhole, FMeshWhole, EV, EF, FE, rawField, matching, effort);

        elapsed = timer.getElapsedTime();
        std::cout << "principal matching elapsed: " << elapsed << std::endl;

        // singularities
        directional::effort_to_indices(VMeshWhole, FMeshWhole, EV, EF, effort, matching, N, singVertices, singIndices);

        elapsed = timer.getElapsedTime();
        std::cout << "singularities elapsed: " << elapsed << std::endl;

        // INTEGRATION

        // integration
        std::cout << "Setting up Integration #" << i << std::endl;

        directional::IntegrationData intData(N);
        directional::setup_integration(VMeshWhole, FMeshWhole, EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);

        elapsed = timer.getElapsedTime();
        std::cout << "Integration setup elapsed: " << elapsed << std::endl;

        intData.verbose = false;
        intData.integralSeamless = true;
        intData.roundSeams = false;

        std::cout << "SOLVING INTEGRATION FOR N=" << N << std::endl;

        directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut, NFunction, NCornerFunction);

        std::cout << "DONE!" << std::endl;

        elapsed = timer.getElapsedTime();
        std::cout << "Integration solve elapsed: " << elapsed << std::endl;

        // MESHING

        // The meshing unit is independent from the integration unit, and can be potentially used with external functions; one should fill the MeshFunctionIsolinesData structure with the input, and call mesh_function_isolines().

        //setting up mesh data from itnegration data
        directional::MeshFunctionIsolinesData mfiData;
        directional::setup_mesh_function_isolines(VMeshCut, FMeshCut, intData, mfiData);

        elapsed = timer.getElapsedTime();
        std::cout << "meshing setup elapsed: " << elapsed << std::endl;

        //meshing
        directional::mesh_function_isolines(VMeshWhole, FMeshWhole, EV, EF, FE, mfiData, verbose, VPolyMesh, DPolyMesh, FPolyMesh);

        elapsed = timer.getElapsedTime();
        std::cout << "meshing elapsed: " << elapsed << std::endl;

        // SAVING

        //saving
        hedra::polygonal_write_OFF(outfilename, VPolyMesh, DPolyMesh, FPolyMesh);

        std::cout << "SAVED!" << std::endl;

        elapsed = timer.getElapsedTime();
        std::cout << "saving elapsed: " << elapsed << std::endl;

        // TIME

        timer.stop();
        elapsed = timer.getElapsedTime();
        std::cout << "Total elapsed: " << elapsed << std::endl;

        // SHOWING

        //raw field mesh
        directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, directional::indexed_glyph_colors(combedField), VField, FField, CField, 1.0);

        if (i == 0)
        {
            viewer.append_mesh();
            viewer.data_list[1].clear();
            viewer.data_list[1].set_mesh(VField, FField);
            viewer.data_list[1].set_colors(CField);
            viewer.data_list[1].show_faces = true;
            viewer.data_list[1].show_lines = false;
        }

        //singularity mesh
        directional::singularity_spheres(VMeshWhole, FMeshWhole, N, singVertices, singIndices, VSings, FSings, CSings, 2.5);

        if (i == 0)
        {
            viewer.append_mesh();
            viewer.data_list[2].clear();
            viewer.data_list[2].set_mesh(VSings, FSings);
            viewer.data_list[2].set_colors(CSings);
            viewer.data_list[2].show_faces = true;
            viewer.data_list[2].show_lines = false;
        }

        //seams mesh
        Eigen::VectorXi isSeam = Eigen::VectorXi::Zero(EV.rows());
        for (int i = 0; i < FE.rows(); i++)
            for (int j = 0; j < 3; j++)
                if (intData.face2cut(i, j))
                    isSeam(FE(i, j)) = 1;
        directional::seam_lines(VMeshWhole, FMeshWhole, EV, combedMatching, VSeams, FSeams, CSeams, 2.5);

        if (i == 0)
        {
            viewer.append_mesh();
            viewer.data_list[3].clear();
            viewer.data_list[3].set_mesh(VSeams, FSeams);
            viewer.data_list[3].set_colors(CSeams);
            viewer.data_list[3].show_faces = true;
            viewer.data_list[3].show_lines = false;
        }

        directional::branched_isolines(VMeshCut, FMeshCut, NFunction, VIso, FIso, CIso);
        if (i == 0)
        {
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
    viewer.data_list[0].show_lines = false;
    update_raw_field_mesh();
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
