#include "Display.hpp"

Display::Display() :
	#ifdef GLFW3
	window(nullptr),
	#endif
	gl2fallback(false)
{
	glfwInit();

	#ifdef GLFW3
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#else // GLFW3
	glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);
	glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
	glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 3);
	glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#endif // GLFW3


	// Open a window and create its OpenGL context
	#ifndef GLFW3
	glfwOpenWindow( 512, 512,6,5,6,0,0,0, GLFW_WINDOW);
	glfwSetWindowTitle("Phase Space View");
	#else // GLFW3
	window = glfwCreateWindow( 512, 512, "Phase Space View", NULL, NULL);
	if( window == nullptr ) {
		glfwTerminate();
		std::cout << "Failed to open OpenGl 3 window, will try OpenGL 2." << std::endl;
		gl2fallback = true;
		glfwInit();
		glfwWindowHint(GLFW_SAMPLES, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
		window = glfwCreateWindow( 512, 512, "Phace Space View", NULL, NULL);

		if( window == nullptr ) {
			std::cerr << "Failed to initialize GLFW" << std::endl;
			glfwTerminate();
		}
	}
	glfwMakeContextCurrent(window);
	#endif // GLFW3

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		std::cerr << "Failed to initialize GLEW" << std::endl;
	}

	// Ensure we can capture the escape key being pressed below
	#ifdef GLFW3
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	#else
	glfwEnable(GLFW_STICKY_KEYS);
	#endif

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	if (!gl2fallback) {
		programID = LoadShaders("gl/gl3.vertexshader","gl/gl3.fragmentshader");
	} else {
		programID = LoadShaders("gl/gl2.vertexshader","gl/gl2.fragmentshader");
	}

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");

	// Projection matrix : 45° Field of View, 1:1 ratio, display range : 0.1 unit <-> 100 units
	glm::mat4 Projection = glm::perspective(45.0f, 1.0f, 0.1f, 100.0f);
	// Camera matrix
	glm::mat4 View	   = glm::lookAt(
								glm::vec3(0.5,0.5,1.0),
								glm::vec3(0.5,0.5,0),
								glm::vec3(0,1,0)
						   );
	// Model matrix : an identity matrix (model will be at the origin)
	glm::mat4 Model	  = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	MVP  = Projection * View * Model;

	static const GLfloat g_vertex_buffer_data[] = {
		0.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		1.0f, 1.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		1.0f, 0.0f, 0.0f
	};

	/* Two UV coordinatesfor each vertex.
	 * Cordinates are switched because Mesh2D has data arranged as [x][y],
	 * which is not the way OpenGL would expect it.
	 */
	static const GLfloat g_uv_buffer_data[] = {
		0.0f, 0.0f,
		0.0f, 1.0f,
		1.0f, 0.0f,
		1.0f, 1.0f,
		1.0f, 0.0f,
		0.0f, 1.0f
	};

	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data),
				 g_vertex_buffer_data, GL_STATIC_DRAW);

	glGenBuffers(1, &uvbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_uv_buffer_data),
				 g_uv_buffer_data, GL_STATIC_DRAW);
}

Display::~Display()
{
	delTexture();

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &uvbuffer);
	glDeleteProgram(programID);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

void Display::createTexture(vfps::PhaseSpace* mesh)
{
	glGenTextures (1, &Texture);
	glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
	glBindTexture(GL_TEXTURE_2D,Texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
	if (std::is_same<vfps::meshdata_t,float>::value) {
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
					 mesh->size(0), mesh->size(1),
					 0, GL_RED,
					 GL_FLOAT, mesh->getData());
	} else {
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
					 mesh->size(0), mesh->size(1),
					 0, GL_RED,
					 GL_INT, mesh->getData());
	}
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,  GL_CLAMP_TO_BORDER );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,  GL_CLAMP_TO_BORDER );
	glBindTexture(GL_TEXTURE_2D, 0);

	// Get a handle for our "myTextureSampler" uniform
	TextureID = glGetUniformLocation(programID, "myTextureSampler");
}

void Display::delTexture()
{
	glDeleteTextures(1, &TextureID);
	glDeleteTextures(1, &Texture);
}

void Display::draw() {
	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Use our shader
	glUseProgram(programID);

	// Send our transformation to the currently bound shader,
	// in the "MVP" uniform
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, Texture);
	// Set our "myTextureSampler" sampler to user Texture Unit 0
	glUniform1i(TextureID, 0);

	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glVertexAttribPointer(
		0,				  // attribute. No particular reason for 0, but must match the layout in the shader.
		3,				  // size
		GL_FLOAT,		   // type
		GL_FALSE,		   // normalized?
		0,				  // stride
		(void*)0			// array buffer offset
	);

	// 2nd attribute buffer : UVs
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glVertexAttribPointer(
		1,								// attribute. No particular reason for 1, but must match the layout in the shader.
		2,								// size : U+V => 2
		GL_FLOAT,						 // type
		GL_FALSE,						 // normalized?
		0,								// stride
		(void*)0						  // array buffer offset
	);
	// Draw the triangle !
	glDrawArrays(GL_TRIANGLES, 0, 2*3); // 12*3 indices starting at 0 -> 12 triangles

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);

	// Swap buffers
	#ifdef GLFW3
	glfwSwapBuffers(window);
	#else
	glfwSwapBuffers();
	#endif
	glfwPollEvents();
}

GLuint Display::LoadShaders(const char* vertex_file_path,
							const char* fragment_file_path){

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if(VertexShaderStream.is_open()){
		std::string Line = "";
		while(getline(VertexShaderStream, Line))
			VertexShaderCode += "\n" + Line;
		VertexShaderStream.close();
	}else{
		std::cout << "Impossible to open" << vertex_file_path
				  << "Are you in the right directory?" << std::endl;
		return 0;
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if(FragmentShaderStream.is_open()){
		std::string Line = "";
		while(getline(FragmentShaderStream, Line))
			FragmentShaderCode += "\n" + Line;
		FragmentShaderStream.close();
	}



	GLint Result = GL_FALSE;
	GLint InfoLogLength;



	// Compile Vertex Shader
	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 1 ){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		std::cout << "Error:" << VertexShaderErrorMessage.data() << std::endl;
	}



	// Compile Fragment Shader
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 1 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		std::cout << "Error:" << FragmentShaderErrorMessage.data() << std::endl;
	}



	// Link the program
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 1 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		std::cout << "Error:" << ProgramErrorMessage.data() << std::endl;
	}

	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

	return ProgramID;
}
