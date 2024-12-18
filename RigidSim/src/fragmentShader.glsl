#version 330 core            // minimal GL version support expected from the GPU

struct LightSource {
  vec3 position;
  vec3 color;
  float intensity;
};
uniform LightSource lightSrc;

struct Material {
  vec3 albedo;
  sampler2D albedoTex;
  int albedoTexLoaded;

  sampler2D normalTex;
  int normalTexLoaded;
};
uniform Material material;

uniform vec3 camPos;

in vec3 fPositionModel;
in vec3 fPosition;
in vec3 fNormal;
in vec2 fTexCoord;

out vec4 colorOut; // shader output: the color response attached to this fragment

uniform mat4 modelMat;
uniform mat3 normMat;

void main() {
  vec3 n = (material.normalTexLoaded == 1) ?
    normalize(normMat*((texture(material.normalTex, fTexCoord).rgb - 0.5)*2.0)) : // colors are in [0,1]^3, and normals are in [-1,1]^3
    normalize(fNormal);

  vec3 radiance = vec3(0, 0, 0);
  vec3 wi = normalize(lightSrc.position - fPosition); // unit vector pointing to the light source
  vec3 Li = lightSrc.color*lightSrc.intensity;
  vec3 albedo = material.albedoTexLoaded==1 ? texture(material.albedoTex, fTexCoord).rgb : material.albedo;

  radiance += Li*albedo*max(dot(n, wi), 0);

  colorOut = vec4(radiance, 1.0); // build an RGBA value from an RGB one
}
