// gravity_2d_bgi.cpp
// 2D approximation of the OpenGL simulation, using graphics.h (WinBGI).
// Left click: spawn an object (hold to "grow" mass with right button).
// Right click while initializing: increase mass.
// WASD / arrow keys: pan camera. +/- to zoom. K to toggle pause. Q to quit.

#include <graphics.h>
#include <conio.h>
#include <windows.h>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>

const double REAL_G = 6.6743e-11; // real G (we will scale)
const double SCALE_G = 6e-6;      // tuned gravitational constant for visible sim
const double TIME_STEP = 0.016;   // ~60 FPS
const double MASS_GROW_RATE = 1.15; // multiply mass when right-click grow
const int WIN_W = 1000;
const int WIN_H = 700;

struct Vec3 {
    double x,y,z;
    Vec3(double X=0,double Y=0,double Z=0):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3&o)const{return Vec3(x+o.x,y+o.y,z+o.z);}
    Vec3 operator-(const Vec3&o)const{return Vec3(x-o.x,y-o.y,z-o.z);}
    Vec3 operator*(double s)const{return Vec3(x*s,y*s,z*s);}
};

struct Object {
    Vec3 pos;
    Vec3 vel;
    double mass;     // kg (sim units)
    double density;  // kg/m^3
    bool initializing=false;
    bool glow=false;
    int color;
    Object(Vec3 p, Vec3 v, double m, double d=5515, int col=WHITE, bool G=false)
        :pos(p),vel(v),mass(m),density(d),glow(G),color(col){}
    double radius_screen(double camScale) const {
        // compute radius in screen pixels using mass->physical volume -> scale down:
        double volume = mass / density;
    double r_phys = cbrt((3.0 * volume) / (4.0 * acos(-1.0))); // meters
        // convert to screen pixels using camScale (higher camScale -> bigger objects)
        double pix = r_phys * 1e-4 * camScale; // tuned factor
        if (pix < 2.0) pix = 2.0;
        return pix;
    }
};

// simple perspective projection to 2D screen
void projectToScreen(const Vec3 &p, const Vec3 &camPos, double camScale, int &sx, int &sy) {
    // camera looks down negative z; use simple perspective: scale/(1 + (z - camZ)/Z0)
    double Z0 = 5000.0; // depth constant controlling perspective strength
    double relZ = (p.z - camPos.z);
    double perspective = camScale / (1.0 + relZ / Z0);
    sx = int((p.x - camPos.x) * perspective + WIN_W/2.0);
    sy = int((p.y - camPos.y) * perspective + WIN_H/2.0);
}

void drawGrid(const std::vector<Object>& objs, const Vec3 &camPos, double camScale) {
    // grid lines in world coords; we draw horizontal lines and vertical lines
    // grid "bend" effect: for each line, we offset y based on masses (approx).
    int divisions = 20;
    double worldSize = 20000.0;
    double step = worldSize / divisions;
    // compute simple mass-influence function (sum of masses / distance)
    for (int xi = -divisions/2; xi <= divisions/2; ++xi) {
        double wx = xi * step;
        // draw vertical lines as series of short segments to allow per-segment offset
        for (double wz = -worldSize/2; wz <= worldSize/2; wz += step/6.0) {
            // compute bending offset (y) due to masses
            double bendY = 0.0;
            for (const auto &o : objs) {
                if (o.initializing) continue;
                double dx = (o.pos.x - wx);
                double dz = (o.pos.z - wz);
                double dist = sqrt(dx*dx + dz*dz) + 1.0;
                // influence ~ mass / dist^2 (scaled)
                bendY += (o.mass * 1e-22) / (dist*dist);
            }
            // convert world points to screen
            int x1, y1;
            double wy = -worldSize * 0.35 + bendY * 1000.0; // tune
            projectToScreen(Vec3(wx, wy, wz), camPos, camScale, x1, y1);
            // next point for small segment
            double wz2 = wz + step/6.0;
            double bendY2 = 0.0;
            for (const auto &o : objs) {
                if (o.initializing) continue;
                double dx = (o.pos.x - wx);
                double dz = (o.pos.z - wz2);
                double dist = sqrt(dx*dx + dz*dz) + 1.0;
                bendY2 += (o.mass * 1e-22) / (dist*dist);
            }
            int x2, y2;
            double wy2 = -worldSize * 0.35 + bendY2 * 1000.0;
            projectToScreen(Vec3(wx, wy2, wz2), camPos, camScale, x2, y2);

            // draw line segment
            setcolor(COLOR(200,200,200));
            line(x1,y1,x2,y2);
        }
    }
    for (int zi = -divisions/2; zi <= divisions/2; ++zi) {
        double wz = zi * step;
        for (double wx = -worldSize/2; wx <= worldSize/2; wx += step/6.0) {
            double bendX = 0.0;
            for (const auto &o : objs) {
                if (o.initializing) continue;
                double dx = (o.pos.x - wx);
                double dz = (o.pos.z - wz);
                double dist = sqrt(dx*dx + dz*dz) + 1.0;
                bendX += (o.mass * 1e-22) / (dist*dist);
            }
            int x1,y1;
            double wy = -worldSize * 0.35 + bendX * 1000.0;
            projectToScreen(Vec3(wx, wy, wz), camPos, camScale, x1, y1);

            double wx2 = wx + step/6.0;
            double bendX2 = 0.0;
            for (const auto &o : objs) {
                if (o.initializing) continue;
                double dx = (o.pos.x - wx2);
                double dz = (o.pos.z - wz);
                double dist = sqrt(dx*dx + dz*dz) + 1.0;
                bendX2 += (o.mass * 1e-22) / (dist*dist);
            }
            int x2,y2;
            double wy2 = -worldSize * 0.35 + bendX2 * 1000.0;
            projectToScreen(Vec3(wx2, wy2, wz), camPos, camScale, x2, y2);

            setcolor(COLOR(200,200,200));
            line(x1,y1,x2,y2);
        }
    }
}

void drawObject(const Object &o, const Vec3 &camPos, double camScale) {
    int sx, sy;
    projectToScreen(o.pos, camPos, camScale, sx, sy);
    double r = o.radius_screen(camScale);
    if (o.glow) {
        // stronger filled circle and halo
        for (int i=4;i>=1;--i) {
            setfillstyle(SOLID_FILL, o.color);
            setcolor(o.color);
            fillellipse(sx, sy, (int)(r*(1.0 + 0.2*i)), (int)(r*(1.0 + 0.2*i)));
        }
        setcolor(o.color);
        setfillstyle(SOLID_FILL, o.color);
        fillellipse(sx, sy, (int)r, (int)r);
    } else {
        setcolor(o.color);
        setfillstyle(SOLID_FILL, o.color);
        fillellipse(sx, sy, (int)r, (int)r);
    }
    // optional label (mass)
    std::ostringstream ss;
    ss << (long long)(o.mass);
    std::string s = ss.str();
    setcolor(WHITE);
    // outtextxy expects a non-const char* in this graphics.h; pass a mutable buffer
    outtextxy(sx + (int)r + 2, sy - 6, &s[0]);
}

int main() {
    int gd = DETECT, gm;
        initgraph(&gd, &gm, nullptr);
    initwindow(WIN_W, WIN_H, nullptr);

    // camera
    Vec3 camPos(0,0, -6000);
    double camScale = 0.8;

    std::vector<Object> objs;
    bool running = true;
    bool pause = false;
    bool leftDown = false, rightDown = false;
    clock_t prevTime = clock();

    // initial objects (example)
    objs.emplace_back(Vec3(-5000, 650, -350), Vec3(0,0,1500), 5.97219e22, 5515, LIGHTCYAN, false);
    objs.emplace_back(Vec3(5000, 650, -350), Vec3(0,0,-1500), 5.97219e22, 5515, LIGHTCYAN, false);
    objs.emplace_back(Vec3(0, 0, -350), Vec3(0,0,0), 1.989e25, 5515, LIGHTRED, true);

    // main loop
    while (running) {
        // time step
        clock_t now = clock();
        double dt = double(now - prevTime) / CLOCKS_PER_SEC;
        if (dt < TIME_STEP) {
            Sleep((DWORD)((TIME_STEP - dt) * 1000));
            now = clock();
            dt = double(now - prevTime) / CLOCKS_PER_SEC;
        }
        prevTime = now;

        // input handling (keyboard)
        if (kbhit()) {
            char ch = getch();
            if (ch == 'q' || ch == 'Q') { running = false; break; }
            if (ch == 'k' || ch == 'K') pause = !pause;
            if (ch == '+') camScale *= 1.1;
            if (ch == '-') camScale /= 1.1;
            if (ch=='w' || ch=='W') camPos.y -= 500;
            if (ch=='s' || ch=='S') camPos.y += 500;
            if (ch=='a' || ch=='A') camPos.x -= 500;
            if (ch=='d' || ch=='D') camPos.x += 500;
        }

        // mouse state
        POINT p;
        GetCursorPos(&p);
        HWND hwnd = GetForegroundWindow();
        ScreenToClient(hwnd, &p);
        int mx = p.x, my = p.y;
        // clamp mouse (some environments may differ)
        mx = std::max(0, std::min(WIN_W-1, mx));
        my = std::max(0, std::min(WIN_H-1, my));

        // left click press
        if (GetAsyncKeyState(VK_LBUTTON) & 0x8000) {
            if (!leftDown) {
                // spawn new object at projected world position (inverse projection)
                // simple inverse: assume z fixed near camera front
                double Zguess = -350.0; // spawn plane
                // convert screen to world (approx)
                double perspective = camScale / (1.0 + (Zguess - camPos.z)/5000.0);
                double wx = (mx - WIN_W/2.0) / perspective + camPos.x;
                double wy = (my - WIN_H/2.0) / perspective + camPos.y;
                objs.emplace_back(Vec3(wx, wy, Zguess), Vec3(0,0,0), 1e22, 5515, WHITE, false);
                objs.back().initializing = true;
                leftDown = true;
            }
        } else {
            if (leftDown) {
                // release: finalize initializing -> launched
                objs.back().initializing = false;
                leftDown = false;
            }
        }

        // right button: while initializing grow last object's mass
        if (GetAsyncKeyState(VK_RBUTTON) & 0x8000) {
            if (!rightDown) rightDown = true;
            if (!objs.empty() && objs.back().initializing) {
                objs.back().mass *= MASS_GROW_RATE;
            }
        } else {
            if (rightDown) rightDown = false;
        }

        // physics update
        if (!pause) {
            // compute pairwise gravity
            for (size_t i = 0; i < objs.size(); ++i) {
                Vec3 acc(0,0,0);
                for (size_t j = 0; j < objs.size(); ++j) {
                    if (i==j) continue;
                    Vec3 diff = objs[j].pos - objs[i].pos;
                    double dist = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z) + 1.0;
                    // scaled gravity: G' * m_j / dist^2
                    double gforce = SCALE_G * objs[j].mass / (dist*dist);
                    acc.x += gforce * (diff.x / dist);
                    acc.y += gforce * (diff.y / dist);
                    acc.z += gforce * (diff.z / dist);
                }
                // integrate velocity
                objs[i].vel.x += acc.x * dt;
                objs[i].vel.y += acc.y * dt;
                objs[i].vel.z += acc.z * dt;
            }

            // update positions
            for (auto &o : objs) {
                o.pos.x += o.vel.x * dt;
                o.pos.y += o.vel.y * dt;
                o.pos.z += o.vel.z * dt;
            }

            // collisions: simple merge if overlap in screen space
            for (size_t i = 0; i < objs.size(); ++i) {
                for (size_t j = i+1; j < objs.size(); ++j) {
                    int sx1, sy1, sx2, sy2;
                    projectToScreen(objs[i].pos, camPos, camScale, sx1, sy1);
                    projectToScreen(objs[j].pos, camPos, camScale, sx2, sy2);
                    double dx = sx1 - sx2;
                    double dy = sy1 - sy2;
                    double d = sqrt(dx*dx + dy*dy);
                    double ri = objs[i].radius_screen(camScale);
                    double rj = objs[j].radius_screen(camScale);
                    if (d < (ri + rj) * 0.8) {
                        // merge smaller into larger (conserve mass & momentum)
                        size_t A = i, B = j;
                        if (objs[A].mass < objs[B].mass) std::swap(A,B);
                        // A has larger mass
                        double totalMass = objs[A].mass + objs[B].mass;
                        objs[A].vel.x = (objs[A].vel.x*objs[A].mass + objs[B].vel.x*objs[B].mass) / totalMass;
                        objs[A].vel.y = (objs[A].vel.y*objs[A].mass + objs[B].vel.y*objs[B].mass) / totalMass;
                        objs[A].vel.z = (objs[A].vel.z*objs[A].mass + objs[B].vel.z*objs[B].mass) / totalMass;
                        objs[A].mass = totalMass;
                        objs.erase(objs.begin() + B);
                        if (B < i) { --i; break; }
                    }
                }
            }
        }

        // Rendering
        cleardevice();
        // Background
        setbkcolor(BLACK);
        // Draw grid
        drawGrid(objs, camPos, camScale);

        // Draw objects
        for (const auto &o : objs) {
            drawObject(o, camPos, camScale);
        }

    // HUD
    setcolor(WHITE);
    std::ostringstream hud;
    hud << "Objects: " << objs.size() << "   Pause: " << (pause ? "YES" : "NO") << "   CamScale: " << camScale;
    // store strings in std::string so we can pass a mutable char* (graphics.h expects char*)
    std::string hudStr = hud.str();
    outtextxy(10,10, &hudStr[0]);
    std::string help1 = "Left click: spawn (hold to reposition), Right click (while initializing): grow mass";
    outtextxy(10, 28, &help1[0]);
    std::string help2 = "WASD: pan  +/-: zoom  K: pause/unpause  Q: quit";
    outtextxy(10, 44, &help2[0]);

        // small sleep to control frame rate
        delay(10);
    }

    closegraph();
    return 0;
}
