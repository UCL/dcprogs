if(NOT WIN32)
  # Newer cmakes can't cope with CTEST_CUSTOM_PRE_TEST having arguments, at
  # least on Mac and Linux. Hence the BS below.
  file(WRITE ${CMAKE_BINARY_DIR}/fakeinstall.cc
       "// Workaround for CTEST_CUSTOM_POST_TEST not allowing arguments \n"
       "#include <stdlib.h> \n"
       "int main() \n"
       "{ \n"
       "  // Because CMake has already expanded DCMAKE_INSTALL_RPATH in cmake_install.cmake\n"
       "  // We have to rerun cmake to reconfigure it and set a different rpath in the fake install\n"
       "  system(\"${CMAKE_COMMAND} .. -DCMAKE_INSTALL_RPATH=${CMAKE_BINARY_DIR}/${TEST_INSTALL_DIRECTORY}/lib\");\n"
       "  system(\"${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${TEST_INSTALL_DIRECTORY} -P cmake_install.cmake\");\n"
       "  return system(\"${CMAKE_COMMAND} .. -UCMAKE_INSTALL_RPATH\");\n"
       "}\n")

  add_executable(fake_test_install ${CMAKE_BINARY_DIR}/fakeinstall.cc)

  file(WRITE ${CMAKE_BINARY_DIR}/CTestCustom.cmake
    "set(CTEST_CUSTOM_PRE_TEST ${CMAKE_BINARY_DIR}/fake_test_install)\n")
else(NOT WIN32)

  file(WRITE ${CMAKE_BINARY_DIR}/CTestCustom.cmake
    "set(CTEST_CUSTOM_PRE_TEST \"\\\"${CMAKE_COMMAND}\\\" "
        "-DCMAKE_INSTALL_PREFIX=${TEST_INSTALL_DIRECTORY} "
        "-P cmake_install.cmake\")\n")
endif(NOT WIN32)
