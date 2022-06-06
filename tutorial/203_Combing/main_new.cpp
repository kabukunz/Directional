#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/combing.h>
#include <directional/directional_viewer.h>

int currF = 0, N;
Eigen::MatrixXi F;
Eigen::MatrixXd V;
Eigen::MatrixXd rawField, combedField, barycenters;
directional::DirectionalViewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::VectorXd effort, combedEffort;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
bool showCombed = false;
bool showSingularities = true;

std::string infilename, rawfieldname, combedfieldname, outfilename;

void update_raw_field_mesh()
{
    Eigen::MatrixXd currField = (showCombed ? combedField : rawField);
    viewer.set_field(currField, directional::DirectionalViewer::indexed_glyph_colors(currField, false));
    viewer.toggle_seams(showCombed);
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer &iglViewer, int key, int modifiers)
{
    switch (key)
    {
        // Select vector
    case '1':
        showCombed = !showCombed;
        update_raw_field_mesh();
        break;

    case 'W':
        if (directional::write_raw_field(combedfieldname, rawField))
            std::cout << "Saved combed raw field" << std::endl;
        else
            std::cout << "Unable to save raw field. Error: " << errno << std::endl;
    }

    return true;
}

int main(int argc, char *argv[])
{
    std::cout << "  1        Toggle raw field/Combed field" << std::endl;
    std::cout << "  W        Save raw field" << std::endl;

    // check args
    if (argc < 2)
    {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " <off-file>" << std::endl;
        return 1;
    }

    N = 4;

    std::string cmdlineFile = argv[1];
    infilename = cmdlineFile + ".off";
    rawfieldname = cmdlineFile + ".rawfield";
    combedfieldname = cmdlineFile + ".combed";
    outfilename = cmdlineFile + "-generated.off";

    igl::readOFF(infilename, V, F);
    directional::read_raw_field(rawfieldname, N, rawField);
    igl::edge_topology(V, F, EV, FE, EF);
    igl::barycenter(V, F, barycenters);

    // computing
    directional::principal_matching(V, F, EV, EF, FE, rawField, matching, effort, singVertices, singIndices);
    directional::combing(V, F, EV, EF, FE, rawField, matching, combedField);
    // to get the (mostly trivial) matching of the combed field
    directional::principal_matching(V, F, EV, EF, FE, combedField, combedMatching, combedEffort, singVertices, singIndices);

    // Mesh setup
    viewer.set_mesh(V, F);
    viewer.toggle_mesh_edges(false);
    update_raw_field_mesh();
    viewer.set_singularities(singVertices, singIndices);
    viewer.set_seams(combedMatching);

    viewer.callback_key_down = &key_down;
    viewer.launch();
}
