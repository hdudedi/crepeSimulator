
// Include VCL library
#include "vcl/vcl.hpp"

// Include common part for exercises
#include "helper_scene.hpp"

// Include exercises
#include "exercises/exercises.hpp"




// ************************************** //
// Global data declaration
// ************************************** //

// Storage for shaders indexed by their names
std::map<std::string,GLuint> shaders;

// General shared elements of the scene such as camera and its controler, visual elements, etc
scene_structure scene;

// The graphical interface. Contains Window object and GUI related variables
gui_structure gui;

// Part specific data - you will specify this object in the corresponding exercise part
scene_exercise exercise;


// ************************************** //
// GLFW event listeners
// ************************************** //

void window_size_callback(GLFWwindow* /*window*/, int width, int height);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_click_callback(GLFWwindow* window, int button, int action, int mods);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);


// ************************************** //
// Start program
// ************************************** //

int main()
{
    // ************************************** //
    // Initialization and data setup
    // ************************************** //

    // Initialize external libraries and window
    initialize_interface(gui);

    // Set GLFW events listener
    glfwSetCursorPosCallback(gui.window, cursor_position_callback );
    glfwSetMouseButtonCallback(gui.window, mouse_click_callback);
    glfwSetWindowSizeCallback(gui.window, window_size_callback);
    glfwSetKeyCallback(gui.window, key_callback);

    load_shaders(shaders);
    setup_scene(scene, gui);


    opengl_debug();
    std::cout<<"*** Setup Data ***"<<std::endl;
    exercise.setup_data(shaders, scene, gui);
    std::cout<<"\t [OK] Data setup"<<std::endl;
    opengl_debug();


    // ************************************** //
    // Animation loop
    // ************************************** //



    std::cout<<"*** Start GLFW loop ***"<<std::endl;
    vcl::glfw_fps_counter fps_counter;
    while( !glfwWindowShouldClose(gui.window) )
    {
        opengl_debug();

        // Clear all color and zbuffer information before drawing on the screen
        clear_screen();opengl_debug();
        // Set a white image texture by default
        glBindTexture(GL_TEXTURE_2D,scene.texture_white);

        // Create the basic gui structure with ImGui
        gui_start_basic_structure(gui,scene, shaders);

        // Perform computation and draw calls for each iteration loop
        exercise.frame_draw(shaders, scene, gui); opengl_debug();


        // Render GUI and update window
        ImGui::End();
        scene.camera_control.update = !(ImGui::IsAnyWindowFocused());
        vcl::imgui_render_frame(gui.window);

        update_fps_title(gui.window, gui.window_title, fps_counter);

        glfwSwapBuffers(gui.window);
        glfwPollEvents();
        opengl_debug();

    }
    std::cout<<"*** Stop GLFW loop ***"<<std::endl;

    // Cleanup ImGui and GLFW
    vcl::imgui_cleanup();

    glfwDestroyWindow(gui.window);
    glfwTerminate();

    return 0;
}

void window_size_callback(GLFWwindow* /*window*/, int width, int height)
{
    glViewport(0, 0, width, height);
    scene.camera.perspective.image_aspect = width / static_cast<float>(height);;
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    scene.camera_control.update_mouse_move(scene.camera, window, float(xpos), float(ypos));
}
void mouse_click_callback(GLFWwindow* window, int button, int action, int mods)
{
    ImGui::SetWindowFocus(nullptr);
    scene.camera_control.update_mouse_click(scene.camera, window, button, action, mods);
}


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS){
        //std::cout<<"touche espace pressee !"<<std::endl;
        exercise.space();
    }

    else if(key == GLFW_KEY_RIGHT && action == GLFW_PRESS ){
        exercise.keypoele(0,scene);
    }
    else if(key == GLFW_KEY_LEFT && action == GLFW_PRESS ){
        exercise.keypoele(1,scene);
    }
    else if(key == GLFW_KEY_DOWN && action == GLFW_PRESS){
        exercise.keypoele(2,scene);
    }
    else if(key == GLFW_KEY_UP && action == GLFW_PRESS ){
        exercise.keypoele(3,scene);
    }
    else if(key == GLFW_KEY_A && action == GLFW_PRESS ){
        exercise.keypoele(4,scene);
    }
    else if(key == GLFW_KEY_D && action == GLFW_PRESS  ){
        exercise.keypoele(5,scene);
    }
    else if(key == GLFW_KEY_W && action == GLFW_PRESS ){
        exercise.keypoele(6,scene);
    }
    else if(key == GLFW_KEY_S && action == GLFW_PRESS  ){
        exercise.keypoele(7,scene);
    }
    else if(key == GLFW_KEY_ENTER && action == GLFW_PRESS  ){
        exercise.begin();
    }

}
