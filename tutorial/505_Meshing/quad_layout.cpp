#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
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

#include <igl/Timer.h>
#include <string>

int N = 4;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams, FIso;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams, VIso;
Eigen::MatrixXd CField, CSeams, CSings, CIso;
Eigen::MatrixXd rawField, combedField;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi DPolyMesh;
Eigen::MatrixXi FPolyMesh;
Eigen::MatrixXd VPolyMesh;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd NFunction, NCornerFunction;
directional::DirectionalViewer viewer;

igl::Timer timer;
std::string infilename, rawfieldname, outfilename;
double elapsed;

typedef enum
{
    FIELD,
    INTEGRATION
} ViewingModes;
ViewingModes viewingMode = FIELD;

void update_viewer()
{
    viewer.toggle_field(false);
    viewer.toggle_singularities(false);
    viewer.toggle_seams(false);
    viewer.toggle_isolines(false);
    
    if (viewingMode == FIELD)
    {
        viewer.toggle_field(true);
        viewer.toggle_singularities(true);
        viewer.toggle_seams(true);
        viewer.toggle_isolines(false);
    }
    else
    {
        viewer.toggle_field(false);
        viewer.toggle_singularities(false);
        viewer.toggle_seams(false);
        viewer.toggle_isolines(true);
    }
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer &viewer, int key, int modifiers)
{
    switch (key)
    {
        // Select vector
    case '1':
        viewingMode = FIELD;
        break;
    case '2':
        viewingMode = INTEGRATION;
        break;
    }
    update_viewer();
    return true;
}


int main(int argc, char *argv[])
{
    // check args
    if (argc < 2)
    {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " <off-file>" << std::endl;
        return 1;
    }

    std::string cmdlineFile = argv[1];
    infilename = cmdlineFile + ".off";
    rawfieldname = cmdlineFile + ".rawfield";
    outfilename = cmdlineFile + "-generated.off";

    igl::readOFF(infilename, VMeshWhole, FMeshWhole);

    // RAW FIELD
    
    // TODO: rawfield generation

    directional::read_raw_field(rawfieldname, N, rawField);
    igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

    std::cout << "  1  Loaded field" << std::endl
              << "  2  Show isoline mesh" << std::endl;

    bool verbose = true;


    // TIMER
    timer.start();

    // COMBING

    // principal matching
    directional::principal_matching(VMeshWhole, FMeshWhole, EV, EF, FE, rawField, matching, effort, singVertices, singIndices);

    elapsed = timer.getElapsedTime();
    std::cout << "principal matching elapsed: " << elapsed << std::endl;

    // INTEGRATION

    // integration
    std::cout << "Setting up Integration" << std::endl;

    directional::IntegrationData intData(N);
    directional::setup_integration(VMeshWhole, FMeshWhole, EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);

    elapsed = timer.getElapsedTime();
    std::cout << "Integration setup elapsed: " << elapsed << std::endl;

    // intData.verbose = false;
    // intData.integralSeamless = true;
    // intData.roundSeams = false;

    // std::cout << "SOLVING INTEGRATION FOR N=" << N << std::endl;

    // directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut, NFunction, NCornerFunction);

    // std::cout << "DONE!" << std::endl;

    // elapsed = timer.getElapsedTime();
    // std::cout << "Integration solve elapsed: " << elapsed << std::endl;

    // // MESHING

    // // The meshing unit is independent from the integration unit, and can be potentially used with external functions; one should fill the MeshFunctionIsolinesData structure with the input, and call mesh_function_isolines().

    // // setting up mesh data from itnegration data
    // directional::MeshFunctionIsolinesData mfiData;
    // directional::setup_mesh_function_isolines(VMeshCut, FMeshCut, intData, mfiData);

    // elapsed = timer.getElapsedTime();
    // std::cout << "meshing setup elapsed: " << elapsed << std::endl;

    // // meshing
    // directional::mesh_function_isolines(VMeshWhole, FMeshWhole, EV, EF, FE, mfiData, verbose, VPolyMesh, DPolyMesh, FPolyMesh);

    // elapsed = timer.getElapsedTime();
    // std::cout << "meshing elapsed: " << elapsed << std::endl;

    // SAVING

    // mesh degree
    Eigen::VectorXi D;
    auto F = FMeshCut;
    D.resize(F.rows());

    for (long i = 0; i < F.rows(); i++)
        D(i) = F.row(i).size();

    Eigen::VectorXi DPolyMesh = D;

    // saving
    hedra::polygonal_write_OFF(outfilename, VMeshCut, DPolyMesh, FMeshCut);

    Eigen::MatrixXd emptyMat;
    // igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-rot-seamless.obj", VMeshCut, FMeshCut, emptyMat,     

    std::cout << "SAVED!" << std::endl;

    elapsed = timer.getElapsedTime();
    std::cout << "saving elapsed: " << elapsed << std::endl;

    // TIME

    timer.stop();
    elapsed = timer.getElapsedTime();
    std::cout << "Total elapsed: " << elapsed << std::endl;

    // SHOWING

    viewer.set_mesh(VMeshWhole, FMeshWhole);
    viewer.set_field(combedField, directional::DirectionalViewer::indexed_glyph_colors(combedField));
    viewer.set_singularities(singVertices, singIndices);
    viewer.set_seams(combedMatching);
    viewer.set_isolines(VMeshCut, FMeshCut, NFunction);

    update_viewer();
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
