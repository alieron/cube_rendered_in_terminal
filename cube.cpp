#include <iostream>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <unistd.h>

// Shortcuts so that I don't have to type std:: all the time
using std::cout;
using std::endl;

// Note: std::max and std::min are left alone since math.h also comes with max and min functions

struct vec3 {
    float x, y, z;

    // Default Constructor
    vec3() { x = y = z = 0; };

    // Three dimensional constructor
    vec3(float xValue, float yValue, float zValue) {
        x = xValue;
        y = yValue;
        z = zValue;
    }

    // Overloaded + operator to add vec3s to vec3
    vec3 operator+(const vec3 &vector) const {
        return vec3(x + vector.x, y + vector.y, z + vector.z);
    }

    // Overloaded - operator to subtract vec3 from vec3
    vec3 operator-(const vec3 &vector) const {
        return vec3(x - vector.x, y - vector.y, z - vector.z);
    }

    // Overloaded % operator to perform the cross product
    vec3 operator%(const vec3 &vector) const {
        return vec3(y * vector.z - z * vector.y, z * vector.x - x * vector.z, x * vector.y - y * vector.x);
    }

    // Overloaded * operator to perfor the dot product
    float operator*(const vec3 &vector) const {
        return x * vector.x + y * vector.y + z * vector.z;
    }

    // Transform the vector such that it has a length of 1 aka becomes a unit vector
    void makeUnit() {
        float scalar = 1 / sqrt(x * x + y * y + z * z);

        x *= scalar;
        y *= scalar;
        z *= scalar;
    }
};

std::ostream &operator<<(std::ostream &stream, const vec3 &vec3) {
    stream << "[ " << vec3.x << ", " << vec3.y << ", " << vec3.z << " ]";
    return stream;
}

struct vec2 {
    float x, y;

    // Default Constructor
    vec2() { x = y = 0; };

    // Three dimensional constructor
    vec2(float xValue, float yValue) {
        x = xValue;
        y = yValue;
    }
};

std::ostream &operator<<(std::ostream &stream, const vec2 &vec2) {
    stream << "[ " << vec2.x << ", " << vec2.y << " ]";
    return stream;
}


/* Setting the canvas size by defining the memory buffers */
const int width = 80;
const int height = 40;

const int char_buffer_size = width * height;

char char_buffer[char_buffer_size];
float z_buffer[char_buffer_size];

/* Commonly used characeters, useful to have as constants */
const char space = ' ';
const char newline = '\n';


struct SquarePlane {
    vec3 vertices[4] = {vertices[0], vertices[1], vertices[2], vertices[3]};
    vec3 normal;

    // Default constructor
    SquarePlane() { vertices[0] = vertices[1] = vertices[2] = vertices[3] = vec3(0,0,0); }

    // Four vertex contructor
    SquarePlane(vec3 &A, vec3 &B, vec3 &C, vec3 &D) {
        vertices[0] = A;
        vertices[1] = B;
        vertices[2] = C;
        vertices[3] = D;

        /* Calculating the unit normal vector for this plane */
        normal = ((vertices[1] - vertices[0]) % (vertices[2] - vertices[0]));
        normal.makeUnit();
    }

    // Directionally sensitive must be (clockwise/anti
    bool edgeFunction(const vec2 &a, const vec2 &b, const vec2 &c) const { 
        return ((c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x) > 0); 
    }    

    // Render/rasterizer, takes the primitives(vertices) and draws the shapes as pixels on the screen(in this case writing to a char buffer)
    void render(const vec3 &light) const {
        // int L = 14 * (normal * light);
        int L = 5 * (normal * light);

        /* Determining the bounding box for each plane */
        int maxX = (width / 2) * (1 + ceil(std::max({vertices[0].x, vertices[1].x, vertices[2].x, vertices[3].x}) / 2));
        int maxY = (height / 2) * (1 + ceil(std::max({vertices[0].y, vertices[1].y, vertices[2].y, vertices[3].y}) / 2));
        
        int minX = (width / 2) * (1 + floor(std::min({vertices[0].x, vertices[1].x, vertices[2].x, vertices[3].x}) / 2));
        int minY = (height / 2) * (1 + floor(std::min({vertices[0].y, vertices[1].y, vertices[2].y, vertices[3].y}) / 2));

        /* Simplified calculation of the z(depth) of each plane for use in the z buffer, same for all points on the plane */
        float maxZ = std::max({vertices[0].z, vertices[1].z, vertices[2].z, vertices[3].z}) + 5;
        
        /* Loop through every pixel in the bounding box */
        for (int py = minY; py < maxY + 1; py++) {
            for (int px = minX; px < maxX + 1; px++) {
                vec2 v0 = vec2((width / 2) * (1 + vertices[0].x / 2), (height / 2) * (1 + vertices[0].y / 2));
                vec2 v1 = vec2((width / 2) * (1 + vertices[1].x / 2), (height / 2) * (1 + vertices[1].y / 2));
                vec2 v2 = vec2((width / 2) * (1 + vertices[2].x / 2), (height / 2) * (1 + vertices[2].y / 2));
                vec2 v3 = vec2((width / 2) * (1 + vertices[3].x / 2), (height / 2) * (1 + vertices[3].y / 2));

                vec2 p = vec2(px, py);

                // Normally this is done on a triangle, I just did the edge function check one more time to use the same algorithm on 
                bool inside = true; 
                inside &= edgeFunction(v0, v1, p); 
                inside &= edgeFunction(v1, v2, p); 
                inside &= edgeFunction(v2, v3, p);
                inside &= edgeFunction(v3, v0, p);

                int index = px + 80 * py;

                if (inside == true && height > py && py > 0 && px > 0 && maxZ > z_buffer[index] && width > px) {
                    z_buffer[index] = maxZ;
                    // char_buffer[index] = ".,`-~:^;=!*#$&"[L > 0 ? L : 0];
                    char_buffer[index] = ".=!*#"[L > 0 ? L : 0];
                }
            }
        }
    }

    void render_debug(const vec3 &light, const float &i) const {
        int L = 5 * (normal * light);

        int maxX = (width / 2) * (1 + ceil(std::max({vertices[0].x, vertices[1].x, vertices[2].x, vertices[3].x}) / 2));
        int maxY = (height / 2) * (1 + ceil(std::max({vertices[0].y, vertices[1].y, vertices[2].y, vertices[3].y}) / 2));
        
        int minX = (width / 2) * (1 + floor(std::min({vertices[0].x, vertices[1].x, vertices[2].x, vertices[3].x}) / 2));
        int minY = (height / 2) * (1 + floor(std::min({vertices[0].y, vertices[1].y, vertices[2].y, vertices[3].y}) / 2));

        float maxZ = std::max({vertices[0].z, vertices[1].z, vertices[2].z, vertices[3].z}) + 5;
        
        for (int py = minY; py < maxY + 1; py++) {
            for (int px = minX; px < maxX + 1; px++) {
                vec2 v0 = vec2((width / 2) * (1 + vertices[0].x / 2), (height / 2) * (1 + vertices[0].y / 2));
                vec2 v1 = vec2((width / 2) * (1 + vertices[1].x / 2), (height / 2) * (1 + vertices[1].y / 2));
                vec2 v2 = vec2((width / 2) * (1 + vertices[2].x / 2), (height / 2) * (1 + vertices[2].y / 2));
                vec2 v3 = vec2((width / 2) * (1 + vertices[3].x / 2), (height / 2) * (1 + vertices[3].y / 2));

                vec2 p = vec2(px, py);

                bool inside = true; 
                inside &= edgeFunction(v0, v1, p); 
                inside &= edgeFunction(v1, v2, p); 
                inside &= edgeFunction(v2, v3, p);
                inside &= edgeFunction(v3, v0, p);

                int index = px + 80 * py;

                if (inside == true && height > py && py > 0 && px > 0 && maxZ > z_buffer[index] && width > px) {
                    z_buffer[index] = maxZ;
                    char_buffer[index] = i+48; // Debug using ASCII numbers to represent the sides
                }
            }
        }

        /* Bresenham Line drawing algorithm, used initially, although it only renders the lines between vertices

        int x0 = (width / 2) * (1 + vertices[0].x / 2);
        int y0 = (height / 2) * (1 + vertices[0].y / 2);

        int x1;
        int y1;

        float D = 1/5;

        // int L = 8 * (normal * light);
        // L = std::max(0, L);


        // For each vertex
        for (int i = 1; i < 5; i++) {
            x1 = (width / 2) * (1 + vertices[i < 4 ? i : 0].x / 2);
            y1 = (height / 2) * (1+ vertices[i < 4 ? i : 0].y / 2);

            // cout << x0 << ", " << y0 << ", " << x1 << ", " << y1 << endl;

            int dx =  abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
            int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1; 
            int err = dx + dy, e2; // error value e_xy
            
            for (;;) {  // loop 
                int index = x0 + 80 * y0;
                if (height > y0 && y0 > 0 && x0 > 0 && width > x0) {
                    // L_buffer[index] = L;
                    // char_buffer[index] = ".,-~:;=!*#$@"[L];
                    char_buffer[index] = '.';
                }
                if (x0 == x1 && y0 == y1) break;
                e2 = 2 * err;
                if (e2 >= dy) { err += dy; x0 += sx; } // e_xy+e_x > 0 
                if (e2 <= dx) { err += dx; y0 += sy; } // e_xy+e_y < 0 
            }

            x0 = x1;
            y0 = y1;
        }*/
    }

    // Relavant rotation transforms using rotation matrices, 
    // pre-calculation of sin(theta) and cos(theta) offer slight optimizations
    void rotateX(float &theta) {
        float sin_theta = sin(theta);
        float cos_theta = cos(theta); 

        for (int i = 0; i < 4; i++) {
            vertices[i] = vec3(vertices[i].x, vertices[i].y * cos_theta - vertices[i].z * sin_theta, vertices[i].y * sin_theta + vertices[i].z * cos_theta);
        }

        normal = vec3(normal.x, normal.y * cos_theta - normal.z * sin_theta, normal.y * sin_theta + normal.z * cos_theta);
    }

    void rotateY(float &theta) {
        float sin_theta = sin(theta);
        float cos_theta = cos(theta); 

        for (int i = 0; i < 4; i++) {
            vertices[i] = vec3(vertices[i].x * cos_theta + vertices[i].z * sin_theta, vertices[i].y, vertices[i].z * cos_theta - vertices[i].x * sin_theta);
        }

        normal = vec3(normal.x * cos_theta + normal.z * sin_theta, normal.y, normal.z * cos_theta - normal.x * sin_theta);
    }

    void rotateZ(float &theta) {
        float sin_theta = sin(theta);
        float cos_theta = cos(theta); 

        for (int i = 0; i < 4; i++) {
            vertices[i] = vec3(vertices[i].x * cos_theta - vertices[i].y * sin_theta, vertices[i].x * sin_theta + vertices[i].y * cos_theta, vertices[i].z);
        }

        normal = vec3(normal.x * cos_theta - normal.y * sin_theta, normal.x * sin_theta + normal.y * cos_theta, normal.z);

    }
};

struct Cube {
    SquarePlane faces[6];

    Cube(SquarePlane &top, SquarePlane &bottom, SquarePlane &front, SquarePlane &back, SquarePlane &right, SquarePlane &left) {
        faces[0] = top;
        faces[1] = bottom;
        faces[2] = front;
        faces[3] = back;
        faces[4] = right;
        faces[5] = left;
    }

    void render(const vec3 &light) const {
        for (int i = 0; i < 6; i++) {
            // faces[i].render_debug(light, i); // Debug mode renders the numbers for the faces
            faces[i].render(light);
        }
        // cout << "\x1b[H";

        // faces[5].render(light); // Render individual faces this way

        for (int k = 0; k < char_buffer_size; k++) {
            cout << (k % width ? char_buffer[k] : newline);
        }
        // cout << endl;
    }

    // Just forwards the rotation of the plane to the vertices, when the verticies rotate, the plane rotates
    void rotateX(float &theta) {
        for (int i = 0; i < 6; i++) {
            faces[i].rotateX(theta);
        }
    }

    void rotateY(float &theta) {
        for (int i = 0; i < 6; i++) {
            faces[i].rotateY(theta);
        }
    }

    void rotateZ(float &theta) {
        for (int i = 0; i < 6; i++) {
            faces[i].rotateZ(theta);
        }
    }
};

/* Deprecated function, moved to within the SquarePlane struct scope

void renderSquarePlane(SquarePlane &plane) {
    float maxX = width / 2 + 30 * plane.maxX;
    float minX = width / 2 + 30 * plane.minX;
    float maxY = height / 2 + 30 * plane.maxY;
    float minY = height / 2 + 30 * plane.minY;
    
    int x0 = (width / 2) * (1 + plane.vertices[0].x / 2);
    int y0 = (height / 2) * (1+ plane.vertices[0].y / 2);
    
    int x1;
    int y1;
    
    for (int i = 1; i < 5; i++) {
        x1 = (width / 2) * (1 + plane.vertices[i < 4 ? i : 0].x / 2);
        y1 = (height / 2) * (1+ plane.vertices[i < 4 ? i : 0].y / 2);
    
        // cout << x0 << ", " << y0 << ", " << x1 << ", " << y1 << endl;

        int dx =  abs (x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = -abs (y1 - y0), sy = y0 < y1 ? 1 : -1; 
        int err = dx + dy, e2; // error value e_xy
       
        for (;;) {  // loop
            if (height > y0 && y0 > 0 && x0 > 0 && width > x0) { char_buffer[x0 + 80 * y0] = '#'; }
            if (x0 == x1 && y0 == y1) break;
            e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; } // e_xy+e_x > 0
            if (e2 <= dx) { err += dx; y0 += sy; } // e_xy+e_y < 0
        }

        x0 = x1;
        y0 = y1;
    }


    for (int k = 0; k < char_buffer_size; k++) {
            cout << (k % width ? char_buffer[k] : newline); 
        }
    cout << endl;
}*/

void createCube() {
    // Vertices A to H
    vec3 A = vec3(1.0, 1.0, 1.0);
    vec3 B = vec3(1.0, 1.0, -1.0);
    vec3 C = vec3(-1.0, 1.0, -1.0);
    vec3 D = vec3(-1.0, 1.0, 1.0);

    vec3 E = vec3(-1.0, -1.0, 1.0);
    vec3 F = vec3(1.0, -1.0, 1.0);
    vec3 G = vec3(1.0, -1.0, -1.0);
    vec3 H = vec3(-1.0, -1.0, -1.0);

    // vec3 A = vec3(30.0, 30.0, 30.0);
    // vec3 B = vec3(30.0, 30.0, -30.0);
    // vec3 C = vec3(-30.0, 30.0, -30.0);
    // vec3 D = vec3(-30.0, 30.0, 30.0);

    // vec3 E = vec3(-30.0, -30.0, 30.0);
    // vec3 F = vec3(30.0, -30.0, 30.0);
    // vec3 G = vec3(30.0, -30.0, -30.0);
    // vec3 H = vec3(-30.0, -30.0, -30.0);

    /* The rasterizer is sensitive to the order of the vertices so some trial and error was conducter */

    // SquarePlane top = SquarePlane(A, B, C, D);
    SquarePlane top = SquarePlane(D, C, B, A); // Reversed

    SquarePlane bottom = SquarePlane(G, H, E, F);
    // SquarePlane bottom = SquarePlane(F, E, H, G); // Reversed

    // SquarePlane front = SquarePlane(B, G, H, C);
    SquarePlane front = SquarePlane(C, H, G, B); // Reversed

    // SquarePlane back = SquarePlane(D, E, F, A);
    SquarePlane back = SquarePlane(A, F, E, D); // Reversed

    // SquarePlane right = SquarePlane(A, F, G, B);
    SquarePlane right = SquarePlane(B, G, F, A); // Reversed

    SquarePlane left = SquarePlane(C, H, E, D);
    // SquarePlane left = SquarePlane(D, E, H, C); // Reversed

    Cube cube = Cube(top, bottom, front, back, right, left);

    // Angle rotated per frame, in radians
    float a = 0.04f;
    float b = 0.02f;

    cout << "\x1b[2J";

    for (;;) {
    // for (int i = 0; i < 1; i++) { // for testing
        memset(char_buffer, space, char_buffer_size);
        memset(z_buffer, 0.0f, char_buffer_size * sizeof(float));
        // memset(L_buffer, 0, char_buffer_size * sizeof(int));

        cout << "\x1b[H";
        
        // Control the light source's location
        cube.render(vec3(0, 1 / sqrt(2), -1 / sqrt(2)));
        // cube.render(vec3(0,0,-1));

        // Control the rotation of the cube
        cube.rotateX(a);
        cube.rotateZ(b);

        // Control the speed of the rotation
        // usleep(5000);
        usleep(50000);
        // usleep(100000);
        // usleep(300000);
    }
}

int main() {
    createCube();
    
    // cout << a << endl;

    // for (int k = 0; k < char_buffer_size; k++) {
    //         cout << (k % width ? char_buffer[k] : newline); 
    //     }
    // cout << endl;
    
}