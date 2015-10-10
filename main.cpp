#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Mesh.h"

int gridX = 600;
int gridY = 600;
int gridZ = 600;

const double fovy = 50.;
const double clipNear = .01;
const double clipFar = 1000.0;
double x = 0;
double y = 0;
double z = -2.5;

std::string path = "/Users/rohansawhney/Desktop/developer/C++/direction-fields/kitten.obj";

Mesh mesh;
double avgEdgeLength = 0;
bool success = true;
bool showBasisCycles = false;
bool setGeneratorIndices = false;
int singularities = 2;
std::vector<double> generatorKs;

void printInstructions()
{
    std::cerr << "space: solve for direction field with new random singularities\n"
              << "k: set generator indices to random values\n"
              << "q: toggle direction field/basis cycles\n"
              << "→/←: increase/decrease singularities\n"
              << "↑/↓: move in/out\n"
              << "w/s: move up/down\n"
              << "a/d: move left/right\n"
              << "escape: exit program\n"
              << std::endl;
}

int randIndex()
{
    return rand() % mesh.vertices.size();
}

void roundNumber(double& r)
{
    r = round(r * 100000.0) / 100000.0;
}

double randNumber()
{
    double r = -M_PI + ((double)rand()/RAND_MAX) * 2 * M_PI;
    roundNumber(r);
    return r;
}

void generateIndexValues(Eigen::VectorXd& ks)
{
    double adjustment = 0;
    double sum = mesh.eulerCharacteristic();
    if (sum == 0) {
        adjustment = randNumber();
        sum = adjustment;
    }
    
    double dx = sum / (double)singularities;
    
    int start = singularities % 2 == 0 ? 0 : 1;
    if (start == 1) ks[0] = dx;
    
    for (int i = start; i < singularities; i += 2) {
        double r = randNumber();
        ks[i] = dx + r;
        ks[i+1] = dx - r;
    }
    
    ks[0] -= adjustment;
}

void assignRandomSingularities()
{
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->k = 0;
    }
    
    if (singularities > 0) {
        Eigen::VectorXd ks(singularities);
        generateIndexValues(ks);
        
        mesh.vertices[randIndex()].k = ks[0];
        int i = 1;
        while (i < singularities) {
            int index = randIndex();
            while (mesh.vertices[index].k != 0) {
                index = randIndex();
            }
            
            mesh.vertices[index].k = ks[i];
            i++;
        }
    }
    
    generatorKs.resize(mesh.generators.size());
    for (size_t i = 0; i < generatorKs.size(); i++) {
        if (setGeneratorIndices) generatorKs[i] = randNumber();
        else generatorKs[i] = 0;
    }
}

void init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);
}

void drawBasisCycles()
{
    // draw primal tree
    glColor4f(0.0, 0.0, 1.0, 0.5);
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
        VertexCIter vp = v->parent;
        glVertex3d(v->position.x(), v->position.y(), v->position.z());
        glVertex3d(vp->position.x(), vp->position.y(), vp->position.z());
    }
    glEnd();
    
    // draw cotree tree
    glColor4f(0.0, 1.0, 0.0, 0.5);
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++){
        Eigen::Vector3d v1 = f->centroid();
        Eigen::Vector3d v2 = f->parent->centroid();
        glVertex3d(v1.x(), v1.y(), v1.z());
        glVertex3d(v2.x(), v2.y(), v2.z());
    }
    glEnd();
    
    // draw generators
    glColor4f(1.0, 0.0, 0.0, 0.5);
    glLineWidth(4.0);
    glBegin(GL_LINES);
    for (int i = 0; i < mesh.generators.size(); i++) {
        for (int j = 0; j < mesh.generators[i].size(); j++) {
            HalfEdgeCIter he = mesh.generators[i][j];
            Eigen::Vector3d v1 = he->face->centroid();
            Eigen::Vector3d v2 = he->flip->face->centroid();
            glVertex3d(v1.x(), v1.y(), v1.z());
            glVertex3d(v2.x(), v2.y(), v2.z());
        }
    }
    
    glEnd();
}

void drawDirectionFields()
{
    // draw direction fields
    glColor4f(0.0, 0.0, 1.0, 0.5);
    for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {

        Eigen::Vector3d start = f->centroid();
        Eigen::Vector3d end = start + f->field*avgEdgeLength * 0.75;
        
        glLineWidth(1.0);
        glBegin(GL_LINES);
        glVertex3d(start.x(), start.y(), start.z());
        glVertex3d(end.x(), end.y(), end.z());
        glEnd();
        
        glPointSize(2.5);
        glBegin(GL_POINTS);
        glVertex3d(end.x(), end.y(), end.z());
        glEnd();
    }
    
    // draw singularities
    glColor4f(1.0, 0.0, 0.0, 0.5);
    glPointSize(4.0);
    glBegin(GL_POINTS);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->k != 0) {
            glVertex3d(v->position.x(), v->position.y(), v->position.z());
        }
    }
    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = (double)viewport[2] / (double)viewport[3];
    gluPerspective(fovy, aspect, clipNear, clipFar);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 0, z, x, y, 0, 0, 1, 0);
    
    if (success) {
        if (showBasisCycles) {
            drawBasisCycles();
            
        } else {
            drawDirectionFields();
        }
    }

    glutSwapBuffers();
}

void keyboard(unsigned char key, int x0, int y0)
{
    switch (key) {
        case 27 :
            exit(0);
        case ' ':
            assignRandomSingularities();
            mesh.solve(generatorKs);
            break;
        case 'a':
            x -= 0.03;
            break;
        case 'd':
            x += 0.03;
            break;
        case 'w':
            y += 0.03;
            break;
        case 's':
            y -= 0.03;
            break;
        case 'k':
            setGeneratorIndices = !setGeneratorIndices;
            break;
        case 'q':
            showBasisCycles = !showBasisCycles;
            break;
    }
    
    glutPostRedisplay();
}

void special(int i, int x0, int y0)
{
    switch (i) {
        case GLUT_KEY_UP:
            z += 0.03;
            break;
        case GLUT_KEY_DOWN:
            z -= 0.03;
            break;
        case GLUT_KEY_LEFT:
            singularities--;
            if (singularities < 0) singularities = 0;
            break;
        case GLUT_KEY_RIGHT:
            singularities++;
            if (singularities > mesh.vertices.size() / 5) singularities = (int)mesh.vertices.size() / 5;
            break;
    }
    
    std::stringstream title;
    title << "Direction Field, Singularities: " << singularities;
    glutSetWindowTitle(title.str().c_str());
    
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    success = mesh.read(path);
    assignRandomSingularities();
    mesh.solve(generatorKs);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        avgEdgeLength += e->length();
    }
    avgEdgeLength /= (double)mesh.edges.size();
    
    printInstructions();
    glutInitWindowSize(gridX, gridY);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    std::stringstream title;
    title << "Direction Field, Singularities: " << singularities;
    glutCreateWindow(title.str().c_str());
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMainLoop();
    
    return 0;
}
