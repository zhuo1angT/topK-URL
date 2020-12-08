cmake_minimum_required(VERSION 3.19.1)

set(CMAKE_CXX_STANDARD 17)

PROJECT(topkurl)

ADD_EXECUTABLE(${PROJECT_NAME} topkurl.cxx)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
