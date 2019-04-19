
#include "cloth.hpp"
#include <random>


#ifdef EXERCISE_CLOTH

using namespace vcl;

int phase=0;
float hauteur=-0.97;
float tempsphase=0;

//      -----------------------------phase1------------------------------------------
// Kernels
float W_poly6(const vec3& pi, const vec3& pj, float H);
vec3 W_poly6_grad(const vec3& pi, const vec3& pj, float H);
vec3 W_spiky_grad(const vec3& pi, const vec3& pj, float H);
int to_Key(const vec3 &p);
int to_Key(int x, int y, int z);

// Convert density to pression using tait equation
float density_to_pression(float rho, float rho0, float stiffness);

// Random value generator
std::default_random_engine generator;
std::uniform_real_distribution<float> distrib(0.0,1.0);
// Counter used to save image on hard drive
int counter_image = 0;

//     -------------------------------------------------------Phase2----------------------
float t=0;

std::vector<vec3> posini;
vcl::mat<3,3> rotpoele;

vec3 spring_force(const vec3& pi, const vec3& pj, float L0, float K, float m)
{
    const vec3 p = pi-pj;
    float L = norm(p);
    const vec3 u = normalize(p);

    const vec3 F = -K/m * (L-L0) * u;
    return F;
}

void collision(std::vector<vcl::vec3>& position, std::vector<vcl::vec3>& vitesse,int k,vcl::vec3 p0,vcl::vec3 n, float r );
void collisionsph(std::vector<vcl::vec3>& position, std::vector<vcl::vec3>& vitesse,int k, vec3 posh, float radius );





//     -----------fonctions communes----------------------------------------------------------------

void scene_exercise::setup_data(std::map<std::string,GLuint>& shaders, scene_structure& , gui_structure& gui)
{
    gui.show_frame_camera = false;

    // Load shaders
    //shaders["open_surface"] = create_shader_program("shaders/mesh_back_illumination/mesh.vert.glsl","shaders/mesh_back_illumination/mesh.frag.glsl");
    shaders["open_surface"] = create_shader_program("shaders/mesh/mesh.vert.glsl","shaders/mesh/mesh.frag.glsl");
    shaders["wireframe_quads"] = create_shader_program("shaders/wireframe_quads/shader.vert.glsl","shaders/wireframe_quads/shader.geom.glsl","shaders/wireframe_quads/shader.frag.glsl");
    shaders["crepe"]=create_shader_program("shaders/mesh_back_illumination/mesh.vert.glsl","shaders/mesh_back_illumination/mesh_crepe.frag.glsl");
    shaders["segment_immediate_mode"] = create_shader_program("shaders/segment_immediate_mode/segment_immediate_mode.vert.glsl","shaders/segment_immediate_mode/segment_immediate_mode.frag.glsl");
    //shaders["mesh"] = create_shader_program("shaders/mesh/mesh.vert.glsl","shaders/mesh/mesh.frag.glsl");

    //Pour phase1
    sphere = mesh_drawable( mesh_primitive_sphere(1.0f));
    text=mesh_drawable(mesh_primitive_quad(vec3(-1,-1,1),vec3(1,-1,1),vec3(-1,1,1)));
    text0=mesh_drawable(mesh_primitive_quad(vec3(-1,-1,1),vec3(1,-1,1),vec3(-1,1,1)));


    //pour transition



    // Load textures
    texture_cloth = texture_gpu(image_load_png("data/cloth/crepe.png"));
    texture_wood = texture_gpu(image_load_png("data/cloth/wood.png"));
    textureText= texture_gpu(image_load_png("data/cloth/text.png"));
    textureText0= texture_gpu(image_load_png("data/cloth/text0.png"));

    // Initialize cloth geometry and particles
    initialize_sph();
    initialize();
    initialize_field_image();

    // Default value for simulation parameters
    user_parameters.K = 25.0f;
    user_parameters.m = 1.0f;
    user_parameters.wind = 10.0f;
    user_parameters.mu = 0.1f;
    //user_parameters.mu = 10.0f;

    // Set collision shapes
    collision_shapes.sphere_p = {0,0.1f,0};
    collision_shapes.sphere_r = 0.2f;
    collision_shapes.ground_height = 0.1f;
    collision_shapes.normalpoele=vec3(0,1,0);

    // Init visual models
    //ground = mesh_drawable(mesh_primitive_quad({-1,hauteur+collision_shapes.ground_height-1e-3f,-1}, {1,collision_shapes.ground_height-1e-3f,-1},{-1,collision_shapes.ground_height-1e-3f,1}));
    ground = mesh_drawable(mesh_primitive_quad({-2,hauteur-0.2f,-2}, {2,hauteur-0.2f,-2},{-2,hauteur-0.2f,2}));
    //ground=mesh_drawable(mesh_primitive_disc(1, vec3(0,collision_shapes.ground_height-1e-3f,0), collision_shapes.normalpoele, 100));
    //ground=mesh_drawable(mesh_primitive_disc(1, vec3(0,hauteur+collision_shapes.ground_height-1e-3f,0), collision_shapes.normalpoele, 100));

    gui_display_texture = true;
    gui_display_wireframe = false;
    gui_param.display_field = true;
    gui_param.display_particles = true;
    gui_param.save_field = false;
}


void scene_exercise::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui)
{
    const float dt = timer.update();
    set_gui();


    tempsphase+=dt;


    if(tempsphase>15.0f && tempsphase<20.0f){
        phase=2;
    }

    else if(tempsphase>20.0f){
        phase=3;
    }


    //std::cout<<tempsphase<<std::endl;
    //std::cout<<t<<std::endl;
    //std::cout<<phase<<std::endl;





     //--------Phase 1--------
    const size_t N_substep = 9;
    if(phase==1){
        // Force constant simulation time step
        float h = dt<=1e-6f? 0.0f :timer.scale*0.0003f;
        for(size_t k_substep=0; k_substep<N_substep; ++k_substep)
        {
            // Update values
            update_density();
            update_pression();
            update_acceleration();


            // Numerical integration
            std::map<int, std::vector<particle_element>> new_map;
            const float damping = 0.5f;
            const size_t N = particles.size();
            for(size_t k=0; k<N; ++k)
            {
                vec3& p = particles[k].p;
                vec3& v = particles[k].v;
                vec3& a = particles[k].a;


                v = (1-h*damping)*v + h*a;
                p = p + h*v;


            }



            // Collision
            const float epsilon = 1e-3f;
            for(size_t k=0; k<N; ++k)
            {
                vec3& p = particles[k].p;
                vec3& v = particles[k].v;


                if( p.y<-1 ) {p.y = -1+epsilon*distrib(generator); v.y *= -0.1f; v.x*=0.5f; v.z*=0.5f;}

                if (p.x * p.x + p.z * p.z > 0.16f) {
                    float cosin = dot(vec2(p.x, p.z)/norm(vec2(p.x, p.z)), vec2(1,0));
                    float cocosin = dot(vec2(p.x, p.z)/norm(vec2(p.x, p.z)), vec2(0,1));
                    p.x = cosin * (0.4f - epsilon*distrib(generator)); p.z = cocosin * (0.4f - epsilon*distrib(generator));
                    v.x *= -0.1f; v.z *= -0.1f;
                }

                new_map[to_Key(p)].push_back(particles[k]);
            }


            old_map = new_map;
        }
        T += 0.1;
    }


    //--------Phase 3
    //le temps de cuisson augmente si la crepe est pres de la poele.
    else if(phase==3){
        // Force constant simulation time step
        float h = dt<=1e-6f? 0.0f :timer.scale*0.1f;
        if(position[0][1]<3*collision_shapes.ground_height){
           t+=dt;
        }

        int N=size_t(std::sqrt(position.size()));

        if( ( !simulation_diverged || force_simulation) && h>0)
        {
            compute_forces();
            numerical_integration(h);

            collision_constraints();                 // Detect and solve collision with other shapes

            //hard_constraints();                      // Enforce hard positional constraints

            normal(position, connectivity, normals); // Update normals of the cloth
            detect_simulation_divergence();          // Check if the simulation seems to diverge
        }


        cloth.data_gpu.update_position(position);
        cloth.data_gpu.update_normal(normals);



        //compute rotationpoele
        vec3 vy=vec3(0,1,0);
        vcl::vec3 v=vcl::cross(vy,collision_shapes.normalpoele);
        float c=dot(vy,collision_shapes.normalpoele);
        rotpoele=rotation_from_axis_angle_mat3(v,acos(c));

    }

    display_elements(shaders, scene, gui);

}


void scene_exercise::numerical_integration(float h)
{
    const size_t NN = position.size();

    for(size_t k=0; k<NN; ++k)
    {
        vec3& p = position[k];
        vec3& v = speed[k];
        const vec3& f = force[k];

        v = v + h*f;
        p = p + h*v;
    }
}



void scene_exercise::display_elements(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    //------------phase 1--------------------------------------
    if(phase==1){
        // Display particles

        if(gui_param.display_particles)
        {
            const size_t N = particles.size();
            sphere.uniform_parameter.scaling = sph_param.h/5.0f;
            sphere.uniform_parameter.color = {1,1,0.9};
            for(size_t k=0; k<N; ++k)
            {
                sphere.uniform_parameter.translation = particles[k].p;
                sphere.draw(shaders["open_surface"],scene.camera);
            }
        }



        // Update field image
        if(gui_param.display_field)
        {
            const size_t im_h = field_image_h.im.height;
            const size_t im_w = field_image_h.im.width;
    //        std::vector<unsigned char>& im_data = field_image.im.data;
            std::vector<unsigned char>& im_data_h = field_image_h.im.data;
    #pragma omp parallel for

            for(size_t ky=0; ky<im_h; ++ky)
            {
                for(size_t kx=0; kx<im_w; ++kx)
                {
                    const float x = 1.0f*kx/(im_w-1.0f)-0.5f;
                    const float y = 0.5f-1.0f*ky/(im_h-1.0f);

                    const float f_h = evaluate_display_field({x,-1.0f,y});
                    const float value_h = 0.5f * std::max(f_h-0.20f,0.0f);

                    float rh = 1;
                    float gh = 1;
                    float bh = 0.7f + 0.3f*std::max(0.0f,1 - value_h);
                    float alpha(255);
    //                if (255 * bh > 254) alpha = 255 - 255*bh;
                    if (x*x +y*y > 0.185f) {
                        alpha = 0;
                        rh = 0.3f;
                        gh = 0.35f;
                        bh = 0.4f;}
    //                else if (255 * bh > 254) {
    //                    rh = 0.3f;
    //                    gh = 0.35f;
    //                    bh = 0.4f;
    //                }
                    else {
    //                    alpha = 1.f - bh * 0.1f;

                        rh +=  - std::max(0.f, T - 30) * 0.002f;
                        gh +=  - std::max(0.f,T - 30) * 0.0024f;
                        bh +=  - std::max(0.f,T - 30) * 0.0025f;

    //                    rh = 0.3f * (1 - alpha) + rh * alpha;
    //                    gh = 0.35f * (1 - alpha) + gh * alpha;
    //                    bh = 0.4f * (1 - alpha) + bh * alpha;
                    }


                    im_data_h[4*(kx+im_w*ky)]   = static_cast<unsigned char>(255*std::max(std::min(rh,1.0f),0.0f));
                    im_data_h[4*(kx+im_w*ky)+1] = static_cast<unsigned char>(255*std::max(std::min(gh,1.0f),0.0f));
                    im_data_h[4*(kx+im_w*ky)+2] = static_cast<unsigned char>(255*std::max(std::min(bh,1.0f),0.0f));

                    im_data_h[4*(kx+im_w*ky)+3] = static_cast<unsigned char>(alpha);
                }
            }




            // Display texture
    //        glBindTexture(GL_TEXTURE_2D, field_image.texture_id);
    //        glTexSubImage2D(GL_TEXTURE_2D, 0, 0,0, GLsizei(im_w), GLsizei(im_h), GL_RGBA, GL_UNSIGNED_BYTE, &im_data[0]);
    //        glGenerateMipmap(GL_TEXTURE_2D);
    //        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    //        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    //        field_image.quad.draw(shaders["mesh"],scene.camera);
            glBindTexture(GL_TEXTURE_2D, field_image_h.texture_id);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0,0, GLsizei(im_w), GLsizei(im_h), GL_RGBA, GL_UNSIGNED_BYTE, &im_data_h[0]);
            glGenerateMipmap(GL_TEXTURE_2D);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            field_image_h.quad.draw(shaders["mesh"],scene.camera);
            glBindTexture(GL_TEXTURE_2D, scene.texture_white);


            // Save texture on hard drive
            if( gui_param.save_field )
            {
                const std::string filename = vcl::zero_fill(std::to_string(counter_image),3);
                image_save_png("output/fluid/file_"+filename+".png",field_image_h.im);
                ++counter_image;
            }
        }

    }

    else if(phase==0){
        //display text
        glBindTexture(GL_TEXTURE_2D, textureText0);
        text0.draw(shaders["open_surface"],scene.camera);
        glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    }

    else if(phase==2){
        //display text
        glBindTexture(GL_TEXTURE_2D, textureText);
        text.draw(shaders["open_surface"],scene.camera);
        glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    }


    //------------phase 3--------------------------------------
    else{
        glEnable( GL_POLYGON_OFFSET_FILL );

        // Display cloth
        if(gui_display_texture)
            glBindTexture(GL_TEXTURE_2D, texture_cloth);
        glPolygonOffset( 1.0, 1.0 );

        //cloth.uniform_parameter.color = color;
        cloth.drawcrepe(shaders["crepe"],scene.camera,t);
        glBindTexture(GL_TEXTURE_2D, scene.texture_white);

        if(gui_display_wireframe)
        {
            glPolygonOffset( 1.0, 1.0 );
            cloth.draw(shaders["wireframe_quads"],scene.camera);
        }





        //display
        //poele.uniform_parameter.scaling=0.03f;
        //poele.draw(shaders["open_surface"],scene.camera);






    }

    // Display ground
    //ground.uniform_parameter.rotation=rotpoele;
    glBindTexture(GL_TEXTURE_2D, texture_wood);
    ground.draw(shaders["open_surface"],scene.camera);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

    pan.uniform_parameter.rotation=rotpoele;
    pan.draw(shaders["mesh"], scene.camera);

}



void scene_exercise::set_gui()
{

    // Can set the speed of the animation
    float scale_min = 0.05f;
    float scale_max = 2.0f;
    ImGui::SliderScalar("Time scale", ImGuiDataType_Float, &timer.scale, &scale_min, &scale_max, "%.2f s");

    ImGui::Checkbox("Display field", &gui_param.display_field);
    ImGui::Checkbox("Display particles", &gui_param.display_particles);
    ImGui::Checkbox("Save field on disk", &gui_param.save_field);

    // Start and stop animation
    if (ImGui::Button("Stop"))
        timer.stop();
    if (ImGui::Button("Start"))
        timer.start();


    // Can set the speed of the animation
    //float scale_min = 0.05f;
    //float scale_max = 2.0f;
    ImGui::SliderScalar("Time scale", ImGuiDataType_Float, &timer.scale, &scale_min, &scale_max, "%.2f s");

    float stiffness_min = 0.1f;
    float stiffness_max = 50.0f;
    ImGui::SliderScalar("Stiffness",ImGuiDataType_Float, &user_parameters.K, &stiffness_min, &stiffness_max, "%.2f s");

    float mu_min = 0.0f;
    float mu_max = 1.0f;
    ImGui::SliderScalar("Damping",ImGuiDataType_Float, &user_parameters.mu, &mu_min, &mu_max, "%.2f s");


    float mass_min = 0.1f;
    float mass_max = 3.0f;
    ImGui::SliderScalar("Mass",ImGuiDataType_Float, &user_parameters.m, &mass_min, &mass_max, "%.2f s");

    float wind_min = 0.0f;
    float wind_max = 100.0f;
    ImGui::SliderScalar("Wind",ImGuiDataType_Float, &user_parameters.wind, &wind_min, &wind_max, "%.2f s");

    ImGui::Checkbox("Wireframe",&gui_display_wireframe);
    ImGui::Checkbox("Texture",&gui_display_texture);


    // Start and stop animation
    if (ImGui::Button("Stop anim"))
        timer.stop();
    if (ImGui::Button("Start anim"))
    {
        if( simulation_diverged )
            force_simulation=true;

        timer.start();
    }

    if (ImGui::Button("Initialize Geometry"))
        initialize();

}


//         ------fonctions phase 1------------------------------------------------------------------------------------------------

// Fill an image with field computed as a distance function to the particles
float scene_exercise::evaluate_display_field(const vcl::vec3& p)
{
    float field = 0.0f;
    const float d = 0.05f*0.6f;
    int key = to_Key(p);
    int step = 2;

    for (int x = -step; x <=step; x++) {
        for (int y = 0; y <= 1; y++) {
            for (int z = -step; z <= step; z++) {
                const std::vector<particle_element> neighbors = old_map[key + to_Key(x, y, z)];
                for (const particle_element& part : neighbors) {
                    {
                        const vec3& pi = part.p;
                        const float r  = norm(p-pi);
                        const float u = r/d;
                        if( u < 4)
                            field += std::exp(-u*u);
                    }
                }
            }
        }
    }
//    const size_t N = particles.size();
//    for(size_t i=0; i<N; ++i)
//    {
//        const vec3& pi = particles[i].p;
//        const float r  = norm(p-pi);
//        const float u = r/d;
//        if( u < 4)
//            field += std::exp(-u*u);
//    }
    return field;
}



// Initialize an image where local density is displayed
void scene_exercise::initialize_field_image()
{
    size_t N = 50; // Image dimension

    field_image_h.quad = mesh_primitive_quad({-0.5,-0.94f,-0.5},{0.5,-0.94f,-0.5},{-0.5,-0.94f,0.5});
    field_image_h.im.width = N;
    field_image_h.im.height = N;
    field_image_h.im.data.resize(4*field_image_h.im.width*field_image_h.im.height);
    field_image_h.texture_id = texture_gpu(field_image_h.im);

    field_image_h.quad.uniform_parameter.shading.ambiant = 1.0f;
    field_image_h.quad.uniform_parameter.shading.diffuse = 0.0f;
    field_image_h.quad.uniform_parameter.shading.specular = 0.0f;
}



float density_to_pression(float rho, float rho0, float stiffness)
{
    const float f = std::max(rho/rho0-1,0.0f);
    return stiffness * std::pow(f,1.7);
}

float W_poly6(const vec3& pi, const vec3& pj, float h)
{
    const float h3 = h*h*h;
    const float alpha = 315/64.0f/3.14f/h3;
    const float r = norm(pi-pj);
    if( r<h )
        return alpha*(1-r/h)*(1-r/h)*(1-r/h);
    else
        return 0.0f;
}

vec3 W_poly6_grad(const vec3& pi, const vec3& pj, float h)
{
    const float h9 = h*h*h*h*h*h*h*h*h;
    const float alpha = 315/64.0f/3.14f/h9;
    const float r = norm(pi-pj);
    if( r<h )
        return -3*2*alpha*(h*h-r*r)*(h*h-r*r)*(pi-pj);
    else
        return {0,0,0};
}

vec3 W_spiky_grad(const vec3& pi, const vec3& pj, float h)
{
    const float h6 = h*h*h*h*h*h;
    const float alpha = 15/3.14f/h6;
    const vec3 u = normalize(pi-pj);
    const float r = norm(pi-pj);
    if( r<h )
        return -3*alpha*(h-r)*(h-r)*u;
    else
        return {0,0,0};
}


int to_Key(const vec3& p) {
    int n = 100;
    int key(static_cast<int>(floorf(p.x * n)) + n * n * static_cast<int>(floorf(p.y * n)) + n*static_cast<int>(floorf(p.z * n)));
    return key;
}

int to_Key(int x, int y, int z) {
    int n = 100;
    return x + n*z + n*n*y;
}

void scene_exercise::initialize_sph()
{
    // Influence distance of a particle (size of the kernel)
    const float h = 0.036f;

    // Rest density (consider 1000 Kg/m^3)
    const float rho0 = 2000.0f;

    // Stiffness (consider ~2000 - used in tait equation)
    const float stiffness = 2000.0f;

    // Viscosity parameter
    const float nu = 20.0f;

    // Total mass of a particle (consider rho0 h^2)
    const float m = rho0*h*h;

    // Initial particle spacing (relative to h)
    const float c = 0.61f;


    // Fill a square with particles
    const float epsilon = 1e-3f;
    int id = 0;
    for(float x=h -0.07f; x<0.07f-h; x=x+c*h)
    {
        for(float y=-0.9f+h; y<-0.0f-h; y=y+c*h)
        {
            for(float z=h -0.07f; z<0.07f-h; z=z+c*h)
            {
            particle_element particle;
            particle.p = {x+epsilon*distrib(generator),y,z+epsilon*distrib(generator)}; // a zero value in z position will lead to a 2D simulation
            particle.id = id;
            id++;
            particles.push_back(particle);
            int key(to_Key(particle.p));
            old_map[key].push_back(particle);
            }
        }
    }


    sph_param.h = h;
    sph_param.rho0 = rho0;
    sph_param.nu = nu;
    sph_param.stiffness = stiffness;
    sph_param.m = m;

    T = 20;

    pan = mesh_load_file_obj("data/pan.obj");
    //pan.uniform_parameter.translation = {-0.0,-0.97,-0.35};
    pan.uniform_parameter.translation = {-0.0,-0.97,-0.70};
    pan.uniform_parameter.color = vec3(0.3f, 0.35f, 0.4f);
    //pan.uniform_parameter.scaling = 0.0021;
    pan.uniform_parameter.scaling = 0.004;
}

void scene_exercise::update_density()
{
    const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
        particles[i].rho = 0.0f;

    for(size_t i=0; i<N; ++i)
    {
        const vec3& pi=particles[i].p;
        int key = to_Key(pi);

        int step = 1;

        for (int x = -step; x <=step; x++) {
            for (int y = -step; y <= step; y++) {
                for (int z = -step; z <= step; z++) {
                    const std::vector<particle_element> neighbors = old_map[key + to_Key(x,y,z)];
                    for (const particle_element& p : neighbors) {
                        const vec3& pj = p.p;
                        const float r = norm (pi - pj);
                        if(r<sph_param.h)
                           particles[i].rho += sph_param.m * W_poly6(pi,pj,sph_param.h);
                    }
                }
            }
        }
    }
}

void scene_exercise::update_pression()
{
    const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
    {
        particles[i].pression = density_to_pression(particles[i].rho,sph_param.rho0, sph_param.stiffness);
    }
}

void scene_exercise::update_acceleration()
{
    // gravity
    const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
    {
        particles[i].a = vec3{0,-9.81f,0};
    }

    // SPH forces
    const float h = sph_param.h;
    const float m = sph_param.m;
#pragma omp parallel for
    for(size_t i=0; i<N; ++i)
    {
        const vec3& pi = particles[i].p;
        int key = to_Key(pi);

        const std::vector<particle_element> neighbors = old_map[key];

        int step = 1;

        for (int x = -step; x <=step; x++) {
            for (int y = -step; y <= step; y++) {
                for (int z = -step; z <= step; z++) {
                    const std::vector<particle_element> neighbors = old_map[key + to_Key(x,y,z)];
                    for (const particle_element& p : neighbors) {
                        if (p.id == particles[i].id)
                            continue;
                        const vec3& pj = p.p;
                        const float r = norm (pi - pj);
                        if(r<h)
                        {
                            const vec3& vi = particles[i].v;
                            const vec3& vj = p.v;

                            const float pression_i = particles[i].pression;
                            const float pression_j = p.pression;

                            const float rho_i = particles[i].rho;
                            const float rho_j = p.rho;

                            const vec3 acc_pression = - m * (pression_i/rho_i/rho_i+pression_j/rho_j/rho_j)*W_poly6_grad(pi,pj,h);
                            const vec3 acc_viscosity = 2*sph_param.nu*m/rho_j * dot(vi-vj,pi-pj)/(norm(pi-pj)+h*h/100.0f)*W_poly6_grad(pi,pj,h);

                            particles[i].a += acc_pression + acc_viscosity;
                        }
                    }
                }
            }
        }
    }

}









//fonctions phase 2----------------------------------------------------------------------------------------------------------

void scene_exercise::hard_constraints()
{
    // Fixed positions of the cloth
    for(const auto& constraints : positional_constraints)
        position[constraints.first] = constraints.second;
}


void scene_exercise::detect_simulation_divergence()
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

    const size_t NN = position.size();
    for(size_t k=0; simulation_diverged==false && k<NN; ++k)
    {
        const float f = norm(force[k]);

        if( f!=f) // detect NaN in force
        {
            std::cout<<"NaN detected in forces"<<std::endl;
            simulation_diverged = true;
        }

        if( f>20.0f ) // detect strong force magnitude
        {
            std::cout<<" **** Warning : Strong force magnitude detected "<<f<<" at vertex "<<k<<" ****"<<std::endl;
            simulation_diverged = true;
        }

        if( position[k].x != position[k].x ) // detect NaN in position
        {
            std::cout<<"NaN detected in positions"<<std::endl;
            simulation_diverged = true;
        }

        if(simulation_diverged==true)
        {
            std::cerr<<" **** Simulation has diverged **** "<<std::endl;
            std::cerr<<" > Stop simulation iterations"<<std::endl;
            timer.stop();
        }
    }
#pragma GCC diagnostic pop
}


void collision(std::vector<vcl::vec3>& position, std::vector<vcl::vec3>& vitesse,int k,vcl::vec3 p0,vcl::vec3 n, float r ){
    float detection = dot(position[k]-p0, n);
    float beta=0.8;
    vec3 v=vitesse[k];
    vec3 vperpen=dot(v,n)*n;
    vec3 vparal=v-vperpen;
    //&& norm(vperpen)>0.05f
    if (detection <= r)
    {
        //position
        float distance=r-detection;
        vec3 pnew=position[k]+distance*n;
        position[k]=pnew;

        //vitesse
        float alpha=0.4;

        //float vparallcare=pow(v.x*v.x+v.y*v.y+v.z*v.z,0.5)-vperpen*vperpen;
        //float vparal=pow(vparallcare,0.5);
        vec3 vnew=alpha*vparal-beta*vperpen;
        vitesse[k]=vnew;

    }

}


// Fill value of force applied on each particle
void scene_exercise::compute_forces()
{
    const size_t NN = force.size();      // Total number of particles
    const int N = size_t(std::sqrt(NN)); // Dimension of the grid in u/v direction

    // Update simulation parameters
    // Adapt parameter to scale with the size of the mesh
    simulation_parameters.m = user_parameters.m / float(NN);
    simulation_parameters.K = user_parameters.K / float(NN);


    // Gravity and damping
    const float m = simulation_parameters.m;   // mass of the particle
    const vec3 g = {0,-9.81f,0};
    for(size_t k=0; k<NN; ++k){
        force[k] = m*g;
        force[k]=force[k]-user_parameters.mu*speed[k];
    }
    //forces related to springs
    /*for(int i=0;i<particles.size()-1;i++){
        particle_element pi=particles[i];
        particle_element pi2=particles[i+1];
        vec3 fi=spring_force(pi2.p,pi.p,1/float(N),K,m);
        forces.push_back(fi);
    }
    const vec3 f_ba = spring_force(pb.p,pa.p,0.5f,K,m); // spring
    const vec3 fg = m*g; // weight
    const vec3 f_cb=spring_force(pc.p,pb.p,0.5f,K,m);
    */
    float k=simulation_parameters.K;
    float l0=simulation_parameters.L0;

    // Note: It can be convenient to navigate through the 2D (u,v) parameters of the grid
    //structural
    for(int ku=0; ku<N; ++ku)
    {
        for(int kv=0; kv<N; ++kv)
        {
            if(ku>0){
                vec3 f2=spring_force(position[kv+N*ku],position[kv+N*(ku-1)],l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f2;
            }
            if(kv>0){
                vec3 f1=spring_force(position[kv+N*ku],position[(kv-1)+N*ku],l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f1;
            }
            if(ku<N-1){
                vec3 f4=spring_force(position[kv+N*ku],position[kv+N*(ku+1)],l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f4;
            }
            if(kv<N-1){
                vec3 f3=spring_force(position[kv+N*ku],position[(kv+1)+N*ku],l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f3;
            }

            // the index in the 1D at parameter (ku,kv) is given by kv+N*ku
            // ex. position[kv+N*ku]
            // ...
        }
    }
    //shearing
    //std::cout<<pow(2,1/2)<<"   "<<pow(2,0.5)<<std::endl;
    for(int ku=0; ku<N; ++ku)
    {
        for(int kv=0; kv<N; ++kv)
        {
            if(ku>0 && kv>0){
                vec3 f1=spring_force(position[kv+N*ku],position[(kv-1)+N*(ku-1)],pow(2,0.5)*l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f1;
            }
            if(kv>0 && ku<N-1){
                vec3 f2=spring_force(position[kv+N*ku],position[(kv-1)+N*(ku+1)],pow(2,0.5)*l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f2;
            }
            if(kv<N-1 && ku<N-1){
                vec3 f3=spring_force(position[kv+N*ku],position[(kv+1)+N*(ku+1)],pow(2,0.5)*l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f3;
            }
            if(kv<N-1 && ku>0){
                vec3 f4=spring_force(position[kv+N*ku],position[(kv+1)+N*(ku-1)],pow(2,0.5)*l0,k,m);
                force[kv+N*ku]=force[kv+N*ku]+f4;
            }

            // the index in the 1D at parameter (ku,kv) is given by kv+N*ku
            // ex. position[kv+N*ku]
            // ...
        }
    }
    //bending
    for(int ku=0; ku<N; ++ku)
    {
        for(int kv=0; kv<N; ++kv)
        {
            if(ku>1){
                vec3 f2=spring_force(position[kv+N*ku],position[kv+N*(ku-2)],2*l0,k/2,m);
                force[kv+N*ku]=force[kv+N*ku]+f2;
            }
            if(kv>1){
                vec3 f1=spring_force(position[kv+N*ku],position[(kv-2)+N*ku],2*l0,k/2,m);
                force[kv+N*ku]=force[kv+N*ku]+f1;
            }
            if(ku<N-2){
                vec3 f4=spring_force(position[kv+N*ku],position[kv+N*(ku+2)],2*l0,k/2,m);
                force[kv+N*ku]=force[kv+N*ku]+f4;
            }
            if(kv<N-2){
                vec3 f3=spring_force(position[kv+N*ku],position[(kv+2)+N*ku],2*l0,k/2,m);
                force[kv+N*ku]=force[kv+N*ku]+f3;
            }
        }
    }
    //3
    for(int ku=0; ku<N; ++ku)
        {
            for(int kv=0; kv<N; ++kv)
            {
                if(ku>2){
                    vec3 f2=spring_force(position[kv+N*ku],position[kv+N*(ku-3)],3*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f2;
                }
                if(kv>2){
                    vec3 f1=spring_force(position[kv+N*ku],position[(kv-3)+N*ku],3*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f1;
                }
                if(ku<N-3){
                    vec3 f4=spring_force(position[kv+N*ku],position[kv+N*(ku+3)],3*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f4;
                }
                if(kv<N-3){
                    vec3 f3=spring_force(position[kv+N*ku],position[(kv+3)+N*ku],3*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f3;
                }
            }
        }

    //4
    for(int ku=0; ku<N; ++ku)
        {
            for(int kv=0; kv<N; ++kv)
            {
                if(ku>3){
                    vec3 f2=spring_force(position[kv+N*ku],position[kv+N*(ku-4)],4*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f2;
                }
                if(kv>3){
                    vec3 f1=spring_force(position[kv+N*ku],position[(kv-4)+N*ku],4*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f1;
                }
                if(ku<N-4){
                    vec3 f4=spring_force(position[kv+N*ku],position[kv+N*(ku+4)],4*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f4;
                }
                if(kv<N-4){
                    vec3 f3=spring_force(position[kv+N*ku],position[(kv+4)+N*ku],4*l0,k/2,m);
                    force[kv+N*ku]=force[kv+N*ku]+f3;
                }


            }
        }


    //shearing2
    /*
        for(int ku=0; ku<N; ++ku)
        {
            for(int kv=0; kv<N; ++kv)
            {
                if(ku>1 && kv>1){
                    vec3 f1=spring_force(position[kv+N*ku],position[(kv-2)+N*(ku-2)],pow(8,0.5)*l0,k,m);
                    force[kv+N*ku]=force[kv+N*ku]+f1;
                }
                if(kv>1 && ku<N-2){
                    vec3 f2=spring_force(position[kv+N*ku],position[(kv-2)+N*(ku+2)],pow(8,0.5)*l0,k,m);
                    force[kv+N*ku]=force[kv+N*ku]+f2;
                }
                if(kv<N-2 && ku<N-2){
                    vec3 f3=spring_force(position[kv+N*ku],position[(kv+2)+N*(ku+2)],pow(8,0.5)*l0,k,m);
                    force[kv+N*ku]=force[kv+N*ku]+f3;
                }
                if(kv<N-2 && ku>1){
                    vec3 f4=spring_force(position[kv+N*ku],position[(kv+2)+N*(ku-2)],pow(8,0.5)*l0,k,m);
                    force[kv+N*ku]=force[kv+N*ku]+f4;
                }

                // the index in the 1D at parameter (ku,kv) is given by kv+N*ku
                // ex. position[kv+N*ku]
                // ...
            }
        }*/



    // You have access to the following parameters
    // const float K = simulation_parameters.K;   // spring stiffness
    // const float L0 = simulation_parameters.L0; // edge length in u and v direction at rest

    //wind
    float wind=user_parameters.wind*0.0015;
    for(size_t k=0; k<NN; ++k){
        vec3 n=normals[k];
        vec3 f=wind*dot(n,vec3(0,0,1))*n;
        force[k]=force[k]+f;
    }

    //rappel au centre de la poele
    float kcenter=5;
    for(size_t k=0; k<NN; ++k){
        //if(true){
        //if(position[k].y>30*collision_shapes.ground_height){
        if(position[k].y<10*collision_shapes.ground_height){
            vec3 posk=posini[k];
            vec3 center=vec3(posk.x,position[k].y,posk.z);
            vec3 f=spring_force(position[k],center,norm(center-position[k]),kcenter,m);
            force[k]=force[k]+f;
        }
    }




}

// Handle detection and response to collision with the shape described in "collision_shapes" variable
void scene_exercise::collision_constraints()
{
    const size_t NN = force.size();
    for(size_t k=0; k<NN; ++k)
    {
        vcl::vec3 p0=vec3(0,hauteur+0.8*collision_shapes.ground_height,0);
        collision(position, speed,k,p0,collision_shapes.normalpoele, collision_shapes.ground_height/10);

    }
#include <random>
    // To fill ...
}

// Initialize the geometrical model
// (this function can be called several times by the user)
void scene_exercise::initialize()
{
    // Number of samples of the model (total number of particles is N_cloth x N_cloth)
    const size_t N_cloth = 16;

    // Rest length (length of an edge)
    simulation_parameters.L0 = 1.0f/float(N_cloth-1);

    // Create cloth mesh in its initial position
    // Horizontal grid of length 1 x 1
    const mesh base_cloth = mesh_primitive_grid(N_cloth,N_cloth,1.0f,1.0f,{-0.5,hauteur,+0.5},{0,1,0});

    //const mesh base_cloth=mesh_primitive_disc(1, {-0.5,1,0}, {0,1,0}, N_cloth);

    // Set particle position from cloth geometry
    position = base_cloth.position;
    posini=base_cloth.position;



    // Set hard positional constraints
    positional_constraints[0] = position[0];
    positional_constraints[N_cloth*(N_cloth-1)] = position[N_cloth*(N_cloth-1)];

    // Init particles data (speed, force)
    speed.resize(0);
    speed.resize(position.size());
    force.resize(0);
    force.resize(position.size());

    //position pour crepe;
    //random
    const size_t NN = force.size();
    for(int i=0;i<NN;i++){
        float rx=0.01f*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float ry=0.01f*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float rz=0.01f*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        position[i]+=vec3(rx,ry,rz);
        //std::cout<<posini[i]<<std::endl;
    }

    for(int i=0;i<NN;i++){
        speed[i]+=vec3(0,0.2,0);
    }


    // Store connectivity and normals
    connectivity = base_cloth.connectivity;
    normals = normal(position,connectivity);



    // Send data to GPU
    cloth.data_gpu.clear();
    cloth = mesh_drawable(base_cloth);
    cloth.uniform_parameter.shading.specular = 0.0f;

    simulation_diverged = false;
    force_simulation = false;

    timer.update();
}






//      ---------------fonctions de commande utilisateur------------------


void scene_exercise::space(){
    //std::cout<<"la aussi ca marche"<<std::endl;
    const size_t NN = force.size();
    for(int i=0;i<NN;i++){
        float rx=0.001f*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float ry=0.001f*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float rz=0.001f*static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        position[i]+=vec3(rx,ry,rz);
    }

    for(int i=0;i<NN;i++){
        speed[i]+=vec3(0,0.25,0);
    }
    //for(int i=NN/2;i<NN;i++){
        //speed[i]+=vec3(0,0.3,0);
    //}
    vec3 vy=vec3(0,1,0);
    vec3 v=vec3(0.2*position[0].y,0.5,0);
    float c=dot(vy,v/norm(v));
    //cloth.uniform_parameter.rotation=rotation_from_axis_angle_mat3(v/norm(v),acos(c));;
}


void scene_exercise::keypoele(int i,scene_structure& scene){
    const size_t NN = force.size();
    vec3 vz=scene.camera.camera_position()-vec3(0,0,0);
    vz=vz-vec3(0,vz.y,0);
    vec3 vx=cross(vz,vec3(0,1,0));
    //std::cout<<vz<<std::endl;
    if(i==0){
        //std::cout<<"right"<<std::endl;
        collision_shapes.normalpoele.x+=0.01;
        collision_shapes.normalpoele=collision_shapes.normalpoele/norm(collision_shapes.normalpoele);
    }
    if(i==1){
        //std::cout<<"left"<<std::endl;
        collision_shapes.normalpoele.x-=0.01;
        collision_shapes.normalpoele=collision_shapes.normalpoele/norm(collision_shapes.normalpoele);
    }
    if(i==2){
        //std::cout<<"down"<<std::endl;
        collision_shapes.normalpoele.z+=0.01;
        collision_shapes.normalpoele=collision_shapes.normalpoele/norm(collision_shapes.normalpoele);
    }
    if(i==3){
        //std::cout<<"up"<<std::endl;
        collision_shapes.normalpoele.z-=0.01;
        collision_shapes.normalpoele=collision_shapes.normalpoele/norm(collision_shapes.normalpoele);
    }
    if(i==5){
        for(int i=0;i<NN;i++){
            //speed[i].x+=0.05;
            speed[i]-=0.03*vx;
        }
    }
    if(i==4){
        for(int i=0;i<NN;i++){
            //speed[i].x-=0.05;
            speed[i]+=0.03*vx;
        }
    }
    if(i==7){
        for(int i=0;i<NN;i++){
            //speed[i].z+=0.05;
            speed[i]+=0.03*vz;
        }

    }
    if(i==6){
        for(int i=0;i<NN;i++){
            //speed[i].z-=0.05;
            speed[i]-=0.03*vz;
        }

    }
}


void scene_exercise::begin(){
    std::cout<<"debut !"<<std::endl;
    tempsphase=0;
    t=0;
    phase=1;
    return;
}






#endif
