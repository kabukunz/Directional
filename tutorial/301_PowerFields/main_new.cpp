#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

Eigen::VectorXi b, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField, representative, bc, barycenters;
Eigen::MatrixXcd powerField;
directional::DirectionalViewer viewer;

int N = 4;
bool normalized = false;
bool zeroPressed = false;

std::string infilename, rawfieldname, outfilename;

void recompute_field()
{
    directional::power_field(V, F, b, bc, N, powerField);
}

void update_raw_field_mesh()
{
    directional::power_to_representative(V, F, powerField, N, representative);
    if (normalized)
        representative.rowwise().normalize();

    directional::representative_to_raw(V, F, representative, N, rawField);
    directional::principal_matching(V, F, EV, EF, FE, rawField, matching, effort, singVertices, singIndices);
    viewer.set_field(rawField);
    viewer.set_singularities(singVertices, singIndices);
}

bool key_up(igl::opengl::glfw::Viewer &viewer, int key, int modifiers)
{
    switch (key)
    {
    case '0':
        zeroPressed = false;
        break;
    }
    return true;
}

bool key_down(igl::opengl::glfw::Viewer &iglViewer, int key, int modifiers)
{
    igl::opengl::glfw::Viewer *iglViewerPointer = &iglViewer;
    directional::DirectionalViewer *directionalViewer = static_cast<directional::DirectionalViewer *>(iglViewerPointer);
    switch (key)
    {
        // Toggle field drawing for easier rotation

    case '0':
        zeroPressed = true;
        break;

        // Reset the constraints
    case 'R':
        b.resize(0);
        bc.resize(0, 3);
        recompute_field();
        update_raw_field_mesh();
        directionalViewer->set_selected_faces(b);
        break;

        // Toggle normalization
    case 'N':
        normalized = !normalized;
        update_raw_field_mesh();
        directionalViewer->set_selected_faces(b);
        break;

    case 'W':
        if (directional::write_raw_field(rawfieldname, rawField))
            std::cout << "Saved raw field" << std::endl;
        else
            std::cout << "Unable to save raw field. Error: " << errno << std::endl;
    }

    return true;
}

// Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer &iglViewer, int key, int modifiers)
{
    igl::opengl::glfw::Viewer *iglViewerPointer = &iglViewer;
    directional::DirectionalViewer *directionalViewer = static_cast<directional::DirectionalViewer *>(iglViewerPointer);
    if (!zeroPressed)
        return false;

    int fid;
    Eigen::Vector3d baryInFace;

    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                 viewer.core().proj, viewer.core().viewport, V, F, fid, baryInFace))
    {

        int i;
        for (i = 0; i < b.rows(); i++)
            if (b(i) == fid)
                break;
        if (i == b.rows())
        {
            b.conservativeResize(b.rows() + 1);
            b(i) = fid;
            bc.conservativeResize(bc.rows() + 1, 3);
        }

        // Compute direction from the center of the face to the mouse
        bc.row(i) = (V.row(F(fid, 0)) * baryInFace(0) +
                     V.row(F(fid, 1)) * baryInFace(1) +
                     V.row(F(fid, 2)) * baryInFace(2) - barycenters.row(fid))
                        .normalized();
        recompute_field();
        update_raw_field_mesh();
        directionalViewer->set_selected_faces(b);
        return true;
    }
    return false;
};

int main(int argc, char *argv[])
{

    std::cout << "  R        Reset the constraints" << std::endl
              << "  N        Toggle field normalization" << std::endl
              << "  W        Save raw field" << std::endl
              << "  0+L-bttn Place constraint pointing from the center of face to the cursor" << std::endl;

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

    // Load mesh
    igl::readOFF(infilename, V, F);
    igl::edge_topology(V, F, EV, FE, EF);
    igl::barycenter(V, F, barycenters);

    b.resize(0);
    bc.resize(0, 3);

    // triangle mesh setup
    viewer.set_mesh(V, F);
    recompute_field();
    update_raw_field_mesh();

    viewer.callback_key_down = &key_down;
    viewer.callback_key_up = &key_up;
    viewer.callback_mouse_down = &mouse_down;
    viewer.launch();
}
