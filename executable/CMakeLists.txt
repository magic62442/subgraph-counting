include_directories(../utility/automorphism)
link_directories(../utility/automorphism)
set (PREPROCESS_SRC
        prepare.cpp)
set (SCOPE_SRC
        scope.cpp)
set (BATCH_SRC
        batch.cpp)
set (FIVE_SRC
        5voc.cpp)
add_executable(preprocess.out ${PREPROCESS_SRC})
add_executable(scope.out ${SCOPE_SRC})
add_executable(batch.out ${BATCH_SRC})
add_executable(5voc.out ${FIVE_SRC})
target_link_libraries(preprocess.out graph counting utility)
target_link_libraries(scope.out graph plan counting utility)
target_link_libraries(batch.out graph plan counting utility)
target_link_libraries(5voc.out graph plan counting utility)
