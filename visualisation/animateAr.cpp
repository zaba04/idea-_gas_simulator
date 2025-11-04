#define GL_SILENCE_DEPRECATION
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <deque>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

double L=3;
int step;
int displayInterval=1;
const double pi=4*atan(1.0);
double radius = 0.12; // ZMNIEJSZONY rozmiar atomów
double minExtent[3], maxExtent[3];
int xWindowSize = 1200, yWindowSize = 800;
GLdouble aspectRatio;
GLdouble fovy, nearClip, farClip;
GLdouble eye[3], center[3], up[3];
GLuint sphereID, configID;
int phi, theta;
int angle=5;
int pos=0;
double Ar[1029];
ifstream datafile;
ifstream outfile;
bool running=false;
bool showInfo=true;

// Parametry układu
int n, N;
double mm, epseps, RR, ff, L_sphere, aa, TT, tautau;
int So, Sd, Sout, Sxyz;

// Dane z out.dat
double currentTime = 0.0;
double currentTemp = 0.0;
double currentEkin = 0.0;
double currentEpot = 0.0;
double currentEtot = 0.0;
double currentPressure = 0.0;

// Historia danych dla wykresów (max 200 punktów)
const int MAX_HISTORY = 200;
deque<double> tempHistory;
deque<double> ekinHistory;
deque<double> epotHistory;
deque<double> etotHistory;
deque<double> timeHistory;

void getParameters(string paramFile){
    ifstream fileOut(paramFile.c_str());
    fileOut>>n>>mm>>epseps>>RR>>ff>>L_sphere>>aa>>TT>>tautau>>So>>Sd>>Sout>>Sxyz;
    fileOut.close();
    N=n*n*n;
}

void readOutData(){
    if(outfile.eof()) return;

    string line;
    if(getline(outfile, line)){
        istringstream iss(line);
        iss >> currentTime >> currentTemp >> currentEkin >> currentEpot >> currentEtot >> currentPressure;

        // Dodaj do historii
        timeHistory.push_back(currentTime);
        tempHistory.push_back(currentTemp);
        ekinHistory.push_back(currentEkin);
        epotHistory.push_back(currentEpot);
        etotHistory.push_back(currentEtot);

    }
}

void makeAtom(GLuint listID, double radius){
    int nTheta=12;
    int nPhi=24;
    glNewList(listID, GL_COMPILE);
    glutSolidSphere(radius, nPhi, nTheta);
    glEndList();
}

void makeCrystal(){
    glNewList(configID, GL_COMPILE);

    glPushMatrix();
    glRotated(phi, 0, 1, 0);
    glPushMatrix();
    glRotated(theta, 1, 0, 0);

    // Rysowanie atomów - JASNY NIEBIESKI
    for (int i=0; i<3*N; i+=3){
        glPushMatrix();
        glTranslated(Ar[i], Ar[i+1], Ar[i+2]);

        // Jasny niebieski z subtelnym gradientem
        float yFactor = (Ar[i+1] + L_sphere) / (2.0 * L_sphere);
        glColor3f(0.4f + 0.2f * yFactor,
                  0.7f + 0.2f * yFactor,
                  1.0f);

        glCallList(sphereID);
        glPopMatrix();
    }

    // NACZYNIE - szklaną kulę z siatką
    glDisable(GL_LIGHTING);

    // Wireframe sfery
    glColor4f(0.2f, 0.5f, 0.8f, 0.4f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutWireSphere(L_sphere, 32, 32);

    // Dodatkowe pionowe i poziome kręgi dla lepszej wizualizacji
    glLineWidth(2.0f);
    glColor4f(0.3f, 0.6f, 0.9f, 0.6f);

    // Równik
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i < 64; i++){
        float angle = 2.0f * pi * i / 64.0f;
        glVertex3f(L_sphere * cos(angle), 0, L_sphere * sin(angle));
    }
    glEnd();

    // Południki
    for(int m = 0; m < 4; m++){
        float baseAngle = m * pi / 2.0f;
        glBegin(GL_LINE_LOOP);
        for(int i = 0; i < 64; i++){
            float angle = 2.0f * pi * i / 64.0f;
            glVertex3f(L_sphere * cos(angle) * cos(baseAngle),
                      L_sphere * sin(angle),
                      L_sphere * cos(angle) * sin(baseAngle));
        }
        glEnd();
    }

    glDisable(GL_BLEND);
    glLineWidth(1.0f);
    glEnable(GL_LIGHTING);

    glPopMatrix();
    glPopMatrix();
    glEndList();
}

void drawText(float x, float y, const char* text, void* font = GLUT_BITMAP_HELVETICA_18){
    glRasterPos2f(x, y);
    for(const char* c = text; *c != '\0'; c++){
        glutBitmapCharacter(font, *c);
    }
}

void drawMiniChart(float x, float y, float width, float height,
                   const deque<double>& data, const char* label,
                   float r, float g, float b, double minVal, double maxVal){
    if(data.empty()) return;

    // Ramka
    glColor3f(0.3f, 0.3f, 0.3f);
    glBegin(GL_LINE_LOOP);
        glVertex2f(x, y);
        glVertex2f(x + width, y);
        glVertex2f(x + width, y + height);
        glVertex2f(x, y + height);
    glEnd();

    // Etykieta
    glColor3f(r, g, b);
    drawText(x + 5, y + height - 15, label, GLUT_BITMAP_HELVETICA_12);

    // Normalizacja danych Y
    double dataMin = minVal;
    double dataMax = maxVal;
    if(dataMin == dataMax){
        dataMin = *min_element(data.begin(), data.end());
        dataMax = *max_element(data.begin(), data.end());
    }
    if(dataMax - dataMin < 0.001) dataMax = dataMin + 0.001;

    // ZMIANA: Rysuj wszystkie punkty proporcjonalnie do całkowitej liczby
    glColor3f(r, g, b);
    glBegin(GL_LINE_STRIP);
    for(size_t i = 0; i < data.size(); i++){
        // X proporcjonalne do pozycji w danych (nie do MAX_HISTORY!)
        float px = x + (width * i) / (data.size() > 1 ? data.size() - 1 : 1);
        float py = y + height * (data[i] - dataMin) / (dataMax - dataMin);
        glVertex2f(px, py);
    }
    glEnd();

    // Znacznik aktualnej pozycji (ostatni punkt)
    if(!data.empty()){
        float lastX = x + width;  // Zawsze na końcu
        float lastY = y + height * (data.back() - dataMin) / (dataMax - dataMin);

        glColor3f(1.0f, 1.0f, 0.0f);
        glPointSize(6.0f);
        glBegin(GL_POINTS);
        glVertex2f(lastX, lastY);
        glEnd();
        glPointSize(1.0f);

        // Pulsująca animacja
        float pulse = 0.5f + 0.5f * sin(currentTime * 10.0f);
        glColor4f(1.0f, 1.0f, 0.0f, 0.3f * pulse);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glPointSize(12.0f);
        glBegin(GL_POINTS);
        glVertex2f(lastX, lastY);
        glEnd();
        glDisable(GL_BLEND);
        glPointSize(1.0f);
    }

    // Wartości min/max
    ostringstream oss;
    glColor3f(0.7f, 0.7f, 0.7f);
    oss.str(""); oss << fixed << setprecision(1) << dataMax;
    drawText(x + width + 5, y + height - 5, oss.str().c_str(), GLUT_BITMAP_HELVETICA_10);
    oss.str(""); oss << fixed << setprecision(1) << dataMin;
    drawText(x + width + 5, y + 5, oss.str().c_str(), GLUT_BITMAP_HELVETICA_10);
}

void drawHUD(){
    if(!showInfo) return;

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, xWindowSize, 0, yWindowSize);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    // Panel informacyjny (lewy górny)
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.0f, 0.0f, 0.0f, 0.8f);
    glBegin(GL_QUADS);
        glVertex2f(10, yWindowSize - 200);
        glVertex2f(280, yWindowSize - 200);
        glVertex2f(280, yWindowSize - 10);
        glVertex2f(10, yWindowSize - 10);
    glEnd();

    // Panel wykresów (prawy)
    glColor4f(0.0f, 0.0f, 0.0f, 0.8f);
    glBegin(GL_QUADS);
        glVertex2f(xWindowSize - 410, yWindowSize - 10);
        glVertex2f(xWindowSize - 10, yWindowSize - 10);
        glVertex2f(xWindowSize - 10, yWindowSize - 510);
        glVertex2f(xWindowSize - 410, yWindowSize - 510);
    glEnd();
    glDisable(GL_BLEND);

    // Tekst informacyjny
    ostringstream oss;

    glColor3f(0.0f, 1.0f, 1.0f);
    drawText(20, yWindowSize - 30, "=== ARGON CRYSTAL MD ===");

    glColor3f(1.0f, 1.0f, 0.0f);
    oss.str(""); oss << "Time: " << fixed << setprecision(4) << currentTime;
    drawText(20, yWindowSize - 55, oss.str().c_str());

    glColor3f(1.0f, 0.5f, 0.0f);
    oss.str(""); oss << "Temp: " << fixed << setprecision(2) << currentTemp << " K";
    drawText(20, yWindowSize - 80, oss.str().c_str());

    glColor3f(0.5f, 1.0f, 0.5f);
    oss.str(""); oss << "E_kin: " << fixed << setprecision(3) << currentEkin;
    drawText(20, yWindowSize - 105, oss.str().c_str());

    glColor3f(0.5f, 0.5f, 1.0f);
    oss.str(""); oss << "E_pot: " << fixed << setprecision(3) << currentEpot;
    drawText(20, yWindowSize - 130, oss.str().c_str());

    glColor3f(1.0f, 0.5f, 1.0f);
    oss.str(""); oss << "E_tot: " << fixed << setprecision(3) << currentEtot;
    drawText(20, yWindowSize - 155, oss.str().c_str());

    glColor3f(0.8f, 0.8f, 0.8f);
    oss.str(""); oss << "Atoms: " << N;
    drawText(20, yWindowSize - 180, oss.str().c_str(), GLUT_BITMAP_HELVETICA_12);

    // WYKRESY
    float chartWidth = 360;
    float chartHeight = 100;
    float chartX = xWindowSize - 400;
    float startY = yWindowSize - 40;

    drawMiniChart(chartX, startY - chartHeight, chartWidth, chartHeight,
                  tempHistory, "Temperature", 1.0f, 0.5f, 0.0f, 0, TT * 2.0);

    drawMiniChart(chartX, startY - 2*chartHeight - 20, chartWidth, chartHeight,
                  ekinHistory, "Kinetic Energy", 0.5f, 1.0f, 0.5f, 0, 0);

    drawMiniChart(chartX, startY - 3*chartHeight - 40, chartWidth, chartHeight,
                  epotHistory, "Potential Energy", 0.5f, 0.5f, 1.0f, 0, 0);

    drawMiniChart(chartX, startY - 4*chartHeight - 60, chartWidth, chartHeight,
                  etotHistory, "Total Energy", 1.0f, 0.5f, 1.0f, 0, 0);

    // Kontrolki
    glColor3f(0.7f, 0.7f, 0.7f);
    drawText(20, 40, "[SPACE] Play/Pause  [H] Toggle HUD  [Arrows] Rotate", GLUT_BITMAP_HELVETICA_12);
    drawText(20, 20, "[+/-] Speed  [R] Reset view  [ESC] Exit", GLUT_BITMAP_HELVETICA_12);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

void takeStep(){
    ++step;
    do{
        for (int i=0; i<3*N; i+=3){
            datafile>>Ar[i]>>Ar[i+1]>>Ar[i+2];
        }
        pos+=N;
    }while (!datafile.eof() && pos%N!=0);

    readOutData();

    if (step%displayInterval==0){
        makeCrystal();
        glutPostRedisplay();
    }
}

void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2], center[0], center[1], center[2], up[0], up[1], up[2]);
    glCallList(configID);
    drawHUD();
    glutSwapBuffers();
}

void reshape(int w, int h){
    xWindowSize = w;
    yWindowSize = h;
    glViewport(0, 0, w, h);
    aspectRatio = w/double(h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovy, aspectRatio, nearClip, farClip);
    glMatrixMode(GL_MODELVIEW);
}

void initView(double *minExtent, double *maxExtent){
    glClearColor(0.05f, 0.05f, 0.15f, 1.0f);

    GLfloat lightDiffuse0[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat lightPosition0[] = {1.0, 1.0, 1.0, 0.0};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition0);

    GLfloat lightDiffuse1[] = {0.5, 0.5, 0.7, 1.0};
    GLfloat lightPosition1[] = {-1.0, -0.5, 1.0, 0.0};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse1);
    glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

    double difExtent[3];
    for (int i=0; i<3; ++i)
        difExtent[i] = maxExtent[i] - minExtent[i];
    double dist = 0;
    for (int i=0; i<3; ++i)
        dist += difExtent[i] * difExtent[i];
    dist = sqrt(dist);

    for (int i=0; i<3; ++i)
        center[i] = minExtent[i] + difExtent[i]/2;
    eye[0] = center[0];
    eye[1] = center[1];
    eye[2] = center[2] + dist;
    up[0] = 0;
    up[1] = 1;
    up[2] = 0;

    nearClip = (dist - difExtent[2]/2)/2;
    farClip = 2*(dist + difExtent[2]/2);
    fovy = difExtent[1]/(dist - difExtent[2]/2)/2;
    fovy = 2*atan(fovy)/pi*180;
    fovy *= 1.2;
}

void keyboard(unsigned char key, int x, int y){
    switch(key){
        case 27:
            exit(0);
            break;
        case ' ':
            if (running){
                glutIdleFunc(NULL);
                running = false;
            } else {
                glutIdleFunc(takeStep);
                running = true;
            }
            break;
        case 'h':
        case 'H':
            showInfo = !showInfo;
            glutPostRedisplay();
            break;
        case 'r':
        case 'R':
            phi = theta = 0;
            glutPostRedisplay();
            break;
        case '+':
        case '=':
            displayInterval = max(1, displayInterval - 1);
            break;
        case '-':
        case '_':
            displayInterval++;
            break;
    }
}

void special(int key, int x, int y){
    switch(key){
        case GLUT_KEY_LEFT:
            phi = (phi - angle) % 360;
            break;
        case GLUT_KEY_RIGHT:
            phi = (phi + angle) % 360;
            break;
        case GLUT_KEY_UP:
            theta = (theta - angle) % 360;
            break;
        case GLUT_KEY_DOWN:
            theta = (theta + angle) % 360;
            break;
    }
    makeCrystal();
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y){
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
        if (running){
            glutIdleFunc(NULL);
            running = false;
        } else {
            glutIdleFunc(takeStep);
            running = true;
        }
    }
}

int main(int argc, char *argv[]){
    if (argc < 4){
        cout << "╔════════════════════════════════════════════════════╗" << endl;
        cout << "║   Enhanced Argon Crystal MD Visualizer v3.0       ║" << endl;
        cout << "╚════════════════════════════════════════════════════╝" << endl;
        cout << "\nUsage: " << argv[0] << " <params> <xyz.dat> <out.dat>\n" << endl;
        cout << "Example: ./animateAr params xyz.dat out.dat\n" << endl;
        return 1;
    }

    getParameters(argv[1]);
    datafile.open(argv[2], ifstream::in);
    outfile.open(argv[3], ifstream::in);

    if(!datafile.is_open()){
        cerr << "Error: Cannot open " << argv[2] << endl;
        return 1;
    }
    if(!outfile.is_open()){
        cerr << "Error: Cannot open " << argv[3] << endl;
        return 1;
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(xWindowSize, yWindowSize);
    glutCreateWindow("Argon Crystal MD - Enhanced with Real-time Charts");

    for (int i=0; i<3; ++i){
        minExtent[i] = -L/2;
        maxExtent[i] = L/2;
    }

    initView(minExtent, maxExtent);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMouseFunc(mouse);

    sphereID = glGenLists(1);
    makeAtom(sphereID, radius);
    configID = glGenLists(1);
    makeCrystal();

    cout << "\n✨ Enhanced Visualization started!" << endl;
    cout << "   Atoms: " << N << " (Light blue spheres)" << endl;
    cout << "   Container: Visible spherical vessel" << endl;
    cout << "   Charts: Real-time parameter tracking" << endl;
    cout << "   Controls: SPACE=play/pause, Arrows=rotate, H=toggle HUD\n" << endl;

    glutMainLoop();

    datafile.close();
    outfile.close();

    return 0;
}