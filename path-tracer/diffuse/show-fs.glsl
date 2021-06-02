#version 300 es 
precision highp float;

in vec2 tex;

out vec4 fragmentColor;

uniform struct {
  sampler2D averagedFrameTexture;
} scene;

void main(void) {
  fragmentColor = texture(scene.averagedFrameTexture, tex);
}