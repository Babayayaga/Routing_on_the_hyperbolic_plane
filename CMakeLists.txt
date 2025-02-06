set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.1...3.23)
project(Routing_on_the_hyperbolic_plane)

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_BUILD_TYPE Release)

if(POLICY CMP0071)
    cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core Qt5)
find_package(LEDA QUIET)

find_package(Qt5 QUIET COMPONENTS Widgets)

if(CGAL_Qt5_FOUND
        AND Qt5_FOUND
        AND (CGAL_Core_FOUND OR LEDA_FOUND))
    # ui files, created with Qt Designer
    qt5_wrap_ui(UIS Main.ui)

    qt5_add_resources(RESOURCE_FILES resources/Main.qrc)

    # cpp files
    add_executable ( Routing_on_the_hyperbolic_plane main.cpp ${CGAL_Qt5_RESOURCE_FILES} ${RESOURCE_FILES} ${UIS}
            include/internal/Qt/Input/PointInput.h
            include/internal/PoincareDiskTraits.h
            include/internal/Qt/Input/PolylineInput.h
            include/internal/Qt/GraphicItems/HyperbolicSegmentGraphicsItem.h
            include/internal/Qt/GraphicItems/RoutingScenarioGraphicsItem.h
            include/internal/RoutingScenario.h
            include/internal/Qt/Input/TriangulationSelectTriangle.h
            include/internal/Qt/Input/QObjectInput.h
            include/internal/Qt/GraphicItems/HyperbolicTriangleGraphicsItem.h
            include/internal/RandomDomainGenerator.h
            include/internal/Triangulation_line_face_iterator.h
            include/internal/Qt/Input/NavigateInput.h
            include/internal/BeltramiKleinTraits.h
    )

    add_executable ( Routing_on_the_hyperbolic_plane_2 main.cpp ${CGAL_Qt5_RESOURCE_FILES} ${RESOURCE_FILES} ${UIS}
            include/internal/Qt/Input/PointInput.h
            include/internal/PoincareDiskTraits.h
            include/internal/Qt/Input/PolylineInput.h
            include/internal/Qt/GraphicItems/HyperbolicSegmentGraphicsItem.h
            include/internal/Qt/GraphicItems/RoutingScenarioGraphicsItem.h
            include/internal/RoutingScenario.h
            include/internal/Qt/Input/TriangulationSelectTriangle.h
            include/internal/Qt/Input/QObjectInput.h
            include/internal/Qt/GraphicItems/HyperbolicTriangleGraphicsItem.h
            include/internal/RandomDomainGenerator.h
            include/internal/Triangulation_line_face_iterator.h
            include/internal/Qt/Input/NavigateInput.h
            include/internal/BeltramiKleinTraits.h
    )

    add_executable ( experiments experiment.cpp
            include/internal/Qt/Input/PointInput.h
            include/internal/PoincareDiskTraits.h
            include/internal/Qt/Input/PolylineInput.h
            include/internal/Qt/GraphicItems/HyperbolicSegmentGraphicsItem.h
            include/internal/Qt/GraphicItems/RoutingScenarioGraphicsItem.h
            include/internal/RoutingScenario.h
            include/internal/Qt/Input/TriangulationSelectTriangle.h
            include/internal/Qt/Input/QObjectInput.h
            include/internal/Qt/GraphicItems/HyperbolicTriangleGraphicsItem.h
            include/internal/RandomDomainGenerator.h
            include/internal/Triangulation_line_face_iterator.h
            include/internal/Qt/Input/NavigateInput.h
            include/internal/BeltramiKleinTraits.h
    )

    target_include_directories(Routing_on_the_hyperbolic_plane PRIVATE ./ ./include)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Routing_on_the_hyperbolic_plane )
    target_link_libraries ( Routing_on_the_hyperbolic_plane CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)
    if(CGAL_Core_FOUND)
        target_link_libraries ( Routing_on_the_hyperbolic_plane CGAL::CGAL_Core)
    else()
        target_link_libraries ( Routing_on_the_hyperbolic_plane ${LEDA_LIBRARIES})
    endif()

    target_include_directories(Routing_on_the_hyperbolic_plane_2 PRIVATE ./ ./include)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Routing_on_the_hyperbolic_plane_2 )
    target_link_libraries ( Routing_on_the_hyperbolic_plane_2 CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)
    if(CGAL_Core_FOUND)
        target_link_libraries ( Routing_on_the_hyperbolic_plane_2 CGAL::CGAL_Core)
    else()
        target_link_libraries ( Routing_on_the_hyperbolic_plane_2 ${LEDA_LIBRARIES})
    endif()

    target_include_directories(experiments PRIVATE ./ ./include)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS experiments )
    target_link_libraries ( experiments CGAL::CGAL CGAL::CGAL_Qt5)

    include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
    cgal_add_compilation_test( Routing_on_the_hyperbolic_plane )

else()
    message("NOTICE: This demo requires CGAL_Core (or LEDA), and Qt5 and will not be compiled.")
endif()
