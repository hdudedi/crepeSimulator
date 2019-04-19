#version 330 core

in struct fragment_data
{
    vec4 position;
    vec4 normal;
    vec4 color;
    vec2 texture_uv;
} fragment;

uniform sampler2D texture_sampler;

out vec4 FragColor;

uniform vec3 camera_position;
uniform vec3 color; // object color
uniform float ambiant  = 0.2;
uniform float diffuse  = 0.8;
uniform float specular = 0.5;
uniform float temps;

vec3 light = camera_position+vec3(+1.0, +1.0, 0.0);

void main()
{
    vec3 n = normalize(fragment.normal.xyz);
    vec3 u = normalize(light-fragment.position.xyz);
    vec3 r = reflect(u,n);
    vec3 t = normalize(fragment.position.xyz-camera_position);

    vec3 colordark=vec3(77.0/255,38.0/255,0.0/255);
    vec3 colorclear=vec3(255.0/255,235.0/255,194.0/255);
    float tfin=30;
    vec3 color=vec3(1,1,1);
    vec3 color2;


    float diffuse_value  = diffuse * abs( dot(u,n) );
    float specular_value = specular * pow( abs( dot(r,t) ), 128.0);

    vec3 white = vec3(1.0);
    float calc=pow(fragment.texture_uv[0]-0.5,2)+pow(fragment.texture_uv[1]-0.5,2);
    if(calc>0.25){
        discard;
    }
    if(temps<tfin){
        color=(colorclear-colordark)*(tfin-temps)/tfin+colordark;
        color2=(colorclear-color)*calc/0.25+color;

    }
    else{
        color=colordark;
        color2=(colorclear-color)*calc/0.25+color;
    }

    vec4 color_texture = texture(texture_sampler, fragment.texture_uv)*vec4(color2,1);
    vec3 c = (ambiant+diffuse_value)*color.rgb*fragment.color.rgb*color_texture.rgb + specular_value*white;

    FragColor = vec4(c, color_texture.a*fragment.color.a);
}


