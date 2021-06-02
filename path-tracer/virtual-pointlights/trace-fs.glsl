#version 300 es
precision highp float;

uniform struct {
  vec3 position;
  mat4 rayDirMatrix;
} camera;

uniform struct {
  mat4 surface;
  mat4 clipper;
  vec4 material;
} quadrics[7];

uniform struct {
  vec4 position;
  vec3 powerDensity;
} lights[2];

uniform struct {
  samplerCube envTexture;
  vec4 randoms[64];
  sampler2D previousFrameTexture;
  float iFrame;
} scene;

in vec2 tex;
in vec4 rayDir;

out vec4 fragmentColor;

const float PI   = 3.14159265358979323846264; // PI
const float PHIG = 1.61803398874989484820459 * 00000.1; // Golden Ratio
const float PIG  = 3.14159265358979323846264 * 00000.1; // PI
const float SQ2G = 1.41421356237309504880169 * 10000.0; // Square Root of Two

float goldRand(in vec3 seed){
    return fract(sin(dot(seed.xy*(seed.z+PHIG), vec2(PHIG, PIG)))*SQ2G);
}

float intersectQuadric(vec4 e, vec4 d, mat4 surface, mat4 clipper){
  float a = dot(d * surface, d);
  float b = dot(d * surface, e) + dot(e * surface, d);
  float c = dot(e * surface, e);

  if(abs(a) < 0.001){
  	float t = - c / b;
  	vec4 h = e + d*t;
  	if(dot(h * clipper, h) > 0.0) {
  	  t = -1.0;
 		}
 		return t;
  }

  float discr = b*b - 4.0*a*c;
  if(discr < 0.0){
    return -1.0;
  }
  float t1 = (-b - sqrt(discr)) / (2.0 * a);
  float t2 = (-b + sqrt(discr)) / (2.0 * a);

  vec4 h1 = e + d * t1;
  if(dot(h1 * clipper, h1) > 0.0) {
    t1 = -1.0;
  }

  vec4 h2 = e + d * t2;
  if(dot(h2 * clipper, h2) > 0.0) {
    t2 = -1.0;
  }

  return (t1<0.0)?t2:((t2<0.0)?t1:min(t1, t2));
}

bool findBestHit(vec4 e, vec4 d, out float bestT, out int bestIndex){
  bestT = 10000.0;
  for(int i=2; i<quadrics.length(); i++){
    float t = intersectQuadric(e, d, quadrics[i].surface, quadrics[i].clipper);
    if(t > 0.0 && t < bestT){
      bestT = t;
      bestIndex = i;
    }
  }
  if(bestT < 9999.0)
    return true;
  else
    return false;
}

vec3 directLighting(vec3 x, vec3 n, vec3 v){
	vec3 reflectedRadiance = vec3(0, 0, 0);
	for(int i=0; i<lights.length(); i++){
    vec3 lightPos = lights[i].position.xyz;
    vec3 lightDiff = lightPos - x * lights[i].position.w;
    float lightDist = length(lightDiff);
    vec3 lightDir = lightDiff / lightDist;//normalize(vec3(1, 1, 1));

    vec4 eShadow = vec4(x + n * 0.01, 1);
    vec4 dShadow = vec4(lightDir, 0);
    float shadowT;
    int shadowIndex;
    if(!findBestHit(eShadow, dShadow, shadowT, shadowIndex) ||
           shadowT * lights[i].position.w > lightDist) {

      vec3 lightPowerDensity = lights[i].powerDensity;
      lightPowerDensity /= lightDist * lightDist;
      vec3 diffuseCoeff = vec3(0.3, 0.3, 0);
      vec3 specularCoeff = vec3(0.0, 0.0, 0.0);
      float shininess = 15.0;
      float cosa = dot(n, lightDir);
      if(cosa < 0.0) {
        cosa = 0.0;
      } else {
        reflectedRadiance += lightPowerDensity * cosa * diffuseCoeff;

        float cosb = dot(n, v);
        vec3 halfway = normalize(v + lightDir);
        float cosd = dot(halfway, n);
        if(cosd < 0.0)
        	cosd = 0.0;
          // lightPowerDensity * cosa * BRDF
        reflectedRadiance += lightPowerDensity * specularCoeff *
        	pow(cosd, shininess) * cosa / max(cosa, cosb);
      }
    }
  }
  return reflectedRadiance;
}

void main(void) {
  vec4 e = vec4(camera.position, 1);
  vec4 d = vec4(normalize(rayDir.xyz), 0);

  vec3 w = vec3(1, 1, 1);

  for(int iBounce=0; iBounce<20; iBounce++){
    float t;
    int i;
    if( findBestHit(e, d, t, i)){
      vec4 hit = e + d * t;

      vec4 gradient = hit * quadrics[i].surface + quadrics[i].surface * hit;
      vec3 normal = normalize(gradient.xyz);
      e = hit;
      e.xyz += normal * 0.01;

      // https://en.wikipedia.org/wiki/Obfuscation_(software)
      // ---------
      // TERULETI HENGER FENYFORRAS 
      // ---------

      // 1. valassz egy pontot a hengeren: 
      vec3 lightSamplePos = vec3(-4,3,0); // kozeppont
      float posNoise = 2.0 * (goldRand(scene.randoms[2 *  iBounce].x * vec3(tex *  d.x *  d.y *  d.z * 1024.0 / 10.0, 1.0)) - 0.5);
      lightSamplePos.z += 2.0 * posNoise; // random pont a henger tengelyen
      
      //random szog a henger tengelyere merolegesen
      vec3 lightSampleDir = vec3(1.0,0.0,0.0);
      
      float randomTheta1 = 6.28318530718 * goldRand(scene.randoms[2 *  iBounce + 1].x * vec3(tex *  d.x *  d.y *  d.z * 1024.0 / 10.0, 1.0));
      lightSampleDir.x = lightSampleDir.x * cos(randomTheta1) - lightSampleDir.y * sin(randomTheta1);
      lightSampleDir.y = lightSampleDir.x * sin(randomTheta1) + lightSampleDir.y * cos(randomTheta1);

      lightSampleDir = normalize(lightSampleDir);
      
      // itt a pont :3
      vec3 samplePoint = lightSamplePos + lightSampleDir;

      // 2. latszik-e az arnyalt pont a fenyforrason mintavetelezett pontbol?
      // -> 'arnyeksugar' kilovese
      vec4 eLight = vec4(samplePoint, 1);
      vec4 dLight = vec4(normalize(e.xyz - samplePoint), 0);
      float lightT;
      int lightIndex;
      if(findBestHit(eLight, dLight, lightT, lightIndex))
      {
        vec4 lightHit = eLight + dLight * lightT;

        float hitDist = length(lightHit.xyz - samplePoint);
        float dist = length(e.xyz - samplePoint);
        
        if(hitDist - 0.01 > dist) 
        {
          
          // mi az a form factor? :c
          float costt = dot(normalize(dLight.xyz), lightSampleDir);
          if(costt < 0.0)
          {
            fragmentColor.rgb += -200.0 * vec3(1,1,0) * (costt / (dist * dist));
          }
        }
      }

      //eredmeny += fenyforras emisszioja * geometriai faktor * hozzajarulas

      if(quadrics[i].material.a > 2.99)     // lightsource
      {
        fragmentColor.rgb = quadrics[i].material.rgb;
      }

      // ez akar az uveg lenni; nem mukodik :c
      if(quadrics[i].material.a > 1.99)     // mirror
      {
        // papir + ceruza modszer
        // vec3 tangentDir = -normalize(d.xyz + normalize(normal)); 
        // float alpha = acos(dot(d.xyz, normal));
        // float beta = asin(sin(alpha)/ 1.5);
        // d.xyz = -normalize(normal) + cos(beta) * tangentDir;

        //https://graphics.cs.wisc.edu/WP/cs559-fall2016/files/2016/12/shirley_chapter_13.pdf
        d.xyz = 
        (d.xyz - normal * dot(d.xyz, normal))/1.5 - 
        normal * sqrt(1.0-(1.0-pow(dot(d.xyz, normal), 2.0))/(1.5*1.5));
      
      }
      else if(quadrics[i].material.a > 0.99)     // mirror
      {
        d.xyz = reflect(d.xyz, normal);
      }
      else                                  // diffuse
      {
        d.xyz = normalize(scene.randoms[iBounce].xyz);
        d.xyz = normalize(normal + d.xyz);
        // a (d.x *  d.y *  d.z)-val valo szorzas eredmenyezte legkellemesebb zajt szamomra
        float perPixelNoise = goldRand(vec3(tex *  d.x *  d.y *  d.z * 1024.0 / 10.0, 1.0)) * 6.28318530718;
        // minden random iranyra egy pixelenkent random forgatas
        d.x = cos(perPixelNoise) * d.x + sin(perPixelNoise) * d.z;
        d.z =-sin(perPixelNoise) * d.x + cos(perPixelNoise) * d.z;
      }
      w *= quadrics[i].material.rgb; //kd * pi
    } else {
      fragmentColor.rgb += w * texture ( scene.envTexture, d.xyz).rgb;
      break;
    }
  }
  fragmentColor =
    texture(scene.previousFrameTexture, tex) * (1.0 - 1.0 / scene.iFrame) +
    fragmentColor * 1.0 / scene.iFrame;

}