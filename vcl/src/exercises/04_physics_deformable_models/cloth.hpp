#pragma once

#include "../../exercises/base_exercise/base_exercise.hpp"

#ifdef EXERCISE_CLOTH

// SPH Particle - Phase 1
struct particle_element
{
    vcl::vec3 p; // Position
    vcl::vec3 v; // Speed
    vcl::vec3 a; // Acceleration

    int id; //for recongnition

    // local density and pression
    float rho;
    float pression;
    particle_element() : p{0,0,0},v{0,0,0},a{0,0,0},rho(0),pression(0) {}
};

// SPH simulation parameters - Phase 2
struct sph_parameters
{
    float h;     // influence distance of a particle
    float rho0;  // rest density
    float m;     // total mass of a particle
    float stiffness; // constant of tait equation (relation density / pression)
    float nu;    // viscosity parameter
};

// Image used to display the water appearance
struct field_display
{
    vcl::image im;           // Image storage on CPU
    GLuint texture_id;       // Texture stored on GPU
    vcl::mesh_drawable quad; // Mesh used to display the texture
};


// User parameters available in the GUI
struct gui_parameters
{
    bool display_field;
    bool display_particles;
    bool save_field;
};





struct user_parameters_structure
{
    float m;    // Global mass (to be divided by the number of particles)
    float K;    // Global stiffness (to be divided by the number of particles)
    float mu;   // Damping
    float wind; // Wind magnitude;
};

struct simulation_parameters_structure
{
    float m;  // mass
    float K;  // spring stiffness
    float L0; // spring rest length
};

// Sphere and ground used for collision
struct collision_shapes_structure
{
    vcl::vec3 sphere_p;
    float sphere_r;
    float ground_height;
    vcl::vec3 normalpoele;
};



struct scene_exercise : base_scene_exercise
{
    //--------------------------Phase 1------------

    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void display(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);

    std::vector<particle_element> particles;
    sph_parameters sph_param;
    float T;

    void update_density();
    void update_pression();
    void update_acceleration();

    float evaluate_display_field(const vcl::vec3& p);

    void initialize_sph();
    void initialize_field_image();
    void set_gui();

    gui_parameters gui_param;
    field_display field_image;
    field_display field_image_h;
    vcl::mesh_drawable sphere;
    vcl::mesh_drawable text;
    vcl::mesh_drawable text0;
    vcl::segments_drawable borders;
    vcl::mesh_drawable pan;
    std::map<int, std::vector<particle_element>> old_map;

    vcl::timer_event timer;






    //--------------------------Phase 2------------
    // Particles parameters
    std::vector<vcl::vec3> position;
    std::vector<vcl::vec3> speed;
    std::vector<vcl::vec3> force;

    // Simulation parameters
    simulation_parameters_structure simulation_parameters; // parameters that user can control directly
    user_parameters_structure user_parameters;             // parameters adjusted with respect to mesh size (not controled directly by the user)

    // Cloth mesh elements
    vcl::mesh_drawable cloth;              // Visual model for the cloth
    std::vector<vcl::vec3> normals;        // Normal of the cloth used for rendering and wind force computation
    std::vector<vcl::index3> connectivity; // Connectivity of the triangular model

    // Parameters of the shape used for collision
    collision_shapes_structure collision_shapes;

    // Store index and position of vertices constrained to have a fixed 3D position
    std::map<int,vcl::vec3> positional_constraints;

    // Textures
    GLuint texture_cloth;
    GLuint texture_wood;
    GLuint textureText;
    GLuint textureText0;

    // Visual elements of the scene
    vcl::mesh_drawable ground;
    vcl::mesh_drawable poele;

    // Gui parameters
    bool gui_display_wireframe;
    bool gui_display_texture;

    // Parameters used to control if the simulation runs when a numerical divergence is detected
    bool simulation_diverged; // Active when divergence is detected
    bool force_simulation;    // Force to run simulation even if divergence is detected




    void initialize();
    void collision_constraints();
    void compute_forces();
    void numerical_integration(float h);
    void detect_simulation_divergence();
    void hard_constraints();



    void display_elements(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);

    void space();
    void begin();
    void keypoele(int i,scene_structure& scene);
};






#endif
