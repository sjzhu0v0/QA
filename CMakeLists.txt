cmake_minimum_required(VERSION 3.10)
project(PairFlow LANGUAGES CXX)

# C++17 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 查找 ROOT（含 Minuit）
find_package(ROOT REQUIRED COMPONENTS Minuit)

# 包含目录
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/macro
    ${CMAKE_CURRENT_SOURCE_DIR}/opt
    ${ROOT_INCLUDE_DIRS}
)

# 启用共享库支持
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- 生成 ROOT 字典 ---
set(DICT_OUT ${CMAKE_CURRENT_BINARY_DIR}/MRootDict.cxx)
add_custom_command(
    OUTPUT ${DICT_OUT}
    COMMAND ${ROOT_rootcling_CMD} -f ${DICT_OUT} -c
            -I${CMAKE_CURRENT_SOURCE_DIR}/include
            -I${CMAKE_CURRENT_SOURCE_DIR}/macro
            ${CMAKE_CURRENT_SOURCE_DIR}/opt/EventData.h
            ${CMAKE_CURRENT_SOURCE_DIR}/opt/LinkDef.h
    DEPENDS opt/EventData.h opt/LinkDef.h
)

# --- 构建共享库 libMRootDict.so ---
add_library(MRootDict SHARED
    opt/EventData.cxx
    ${DICT_OUT}
)
target_link_libraries(MRootDict ${ROOT_LIBRARIES})

# --- 所有可执行文件列表 ---
set(SIMPLE_EXES
    process/track/TrackInfoMC.cpp
    macro/event/MultRaw.cpp
    macro/event/MultCalib.cpp
    macro/event/MultPileupCut.cpp
    macro/event/MultQA_AllCut.cpp
    macro/event/MultREFRaw.cpp
)

set(DICT_EXES
    macro/event_jpsi/NJpsiCandidatePerEvent.cpp
    macro/event_jpsi/JpsiAsso.cpp
    macro/event_jpsi/JpsiAsso_BS.cpp
    macro/jpsi/JpsiQA.cpp
    macro/event_jpsi/EventMixingJpsiAsso.cpp
    macro/event_jpsi/MixEventReading.cpp
)

# --- 添加普通可执行文件 ---
foreach(src ${SIMPLE_EXES})
    get_filename_component(target ${src} NAME_WE)
    add_executable(${target}.exe ${src})
    target_link_libraries(${target}.exe ${ROOT_LIBRARIES} Minuit)
endforeach()

# --- 添加依赖 MRootDict 的可执行文件 ---
foreach(src ${DICT_EXES})
    get_filename_component(target ${src} NAME_WE)
    add_executable(${target}.exe ${src})
    target_link_libraries(${target}.exe MRootDict ${ROOT_LIBRARIES} Minuit)
    target_include_directories(${target}.exe PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# --- 统一输出到 bin/ ---
set_target_properties(
    TrackInfoMC.exe MultRaw.exe MultCalib.exe MultPileupCut.exe
    MultQA_AllCut.exe MultREFRaw.exe NJpsiCandidatePerEvent.exe
    JpsiAsso.exe JpsiAsso_BS.exe JpsiQA.exe
    EventMixingJpsiAsso.exe MixEventReading.exe
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# --- Linux: 设置 RPATH 以便运行时找到 libMRootDict.so ---
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set_target_properties(MRootDict PROPERTIES INSTALL_RPATH "$ORIGIN")
    foreach(tgt TrackInfoMC.exe MultRaw.exe MultCalib.exe MultPileupCut.exe
                MultQA_AllCut.exe MultREFRaw.exe NJpsiCandidatePerEvent.exe
                JpsiAsso.exe JpsiAsso_BS.exe JpsiQA.exe
                EventMixingJpsiAsso.exe MixEventReading.exe)
        set_target_properties(${tgt} PROPERTIES INSTALL_RPATH "$ORIGIN")
    endforeach()
endif()
