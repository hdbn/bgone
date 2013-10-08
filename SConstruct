#
# This file is a scons build script for lz factorization programs
# BGS, iBGS, BGL, iBGL, BGT, iBGT
# Copyright 2012 Hideo Bannai & Keisuke Goto
# 

SRC_DIR = "src/"
OUT_DIR = "out/"

import os, sys, glob
import types
env = Environment(ENV = {'PATH' : os.environ['PATH']},
                  CC="gcc",CXX="g++",
                  CCFLAGS="-Wall -fast -DNDEBUG -msse4.2",
                  LINKFLAGS="-Ofast -msse4.2")
# env = Environment(CC="gcc", CXX="/opt/local/bin/g++", 
#                   CCFLAGS="-Wall -fast -DNDEBUG -msse4.2",
#                   LINKFLAGS="-fast -msse4.2"
#                   )
# env = Environment(CC="gcc", CXX="g++", 
#                   CCFLAGS="-Wall -fast -DNDEBUG -msse4.2",
#                   LINKFLAGS="-fast -msse4.2"
#                   )
envDebug = Environment(ENV = {'PATH' : os.environ['PATH']},
    CC="gcc", CXX="g++",
    CCFLAGS="-pg -g -Wall"
    )


def makeprog(onlyfor_sources, progname):
    name, ext = os.path.splitext(os.path.basename(onlyfor_sources))
    # print name, ext
    objects_onlyfor = makeObj(env, onlyfor_sources, '.o')
    env.Program(OUT_DIR + progname, objects_common + objects_onlyfor)


def makeObj(env, srcs, suf):
    if type(srcs) != types.ListType:
        srcs = [srcs]
    ret = []
    for src in srcs:
        name, ext = os.path.splitext(os.path.basename(src))
        # print name, ext
        ret.append(env.Object(OUT_DIR + name+suf, src))
    return ret
# def makeDebugObj(env, srcs):
#     ret = []
#     for src in srcs:
#         name, ext = os.path.splitext(os.path.basename(src))
#         # print name, ext
#         ret.append(env.Object(OUT_DIR + name+'.debug.o', src))
#     return ret


sources_common = [SRC_DIR + src for src in ['bgCommon.cpp',
                                            'divsufsort.c',
                                            'sa2phi.cpp',
                                            'saca-k.cpp']]
objects_common = makeObj(env, sources_common, '.o')

debug_obj_common = makeObj(envDebug, sources_common, '.debug.o')

def makeprogDebug(env, onlyfor_sources, progname, srcs):
    obj_onlyfor = makeObj(env, onlyfor_sources, '.debug.o')
    return env.Program(OUT_DIR + progname + '.debug', obj_onlyfor + srcs)

progs = [
    ['sa2phi_test.cpp', 'sa2phi_test'],
    ['BGtwoMain.cpp', 'BGtwo'],
    ['BGoneMain.cpp', 'BGone']
    # ['sa2phi_test.cpp', 'sa2phi_test'],
    # ['bgt9Main.cpp', 'lzBGT9'],
    # ['lzShuffleMain.cpp', 'lzShuffle']

    ]

for fin, fout in progs:
    # print 'fin',fin
    makeprog(SRC_DIR + fin, fout)
    makeprogDebug(envDebug, [SRC_DIR + fin], fout, debug_obj_common)

# tests: uses google-test
####################################################
#
def runUnitTest(env,target,source):
    import subprocess
    app = str(source[0].abspath)
    if not subprocess.call(app,cwd="./tests"):
        open(str(target[0]),'w').write("PASSED\n")

if os.path.exists('/opt/local/lib/libgtest.a') and os.path.exists('/opt/local/lib/libgtest_main.a'):
    testProg = OUT_DIR + 'runTests'
    envTEST = env.Clone(CPPPATH = ['./', '/opt/local/include'],
                        LIBPATH=['/opt/local/lib', './'],
                        LIBS=['gtest', 'gtest_main', 
                              ])
    envTESTDebug = envDebug.Clone(CPPPATH = ['./', '/opt/local/include'],
                        LIBPATH=['/opt/local/lib', './'],
                        LIBS=['gtest', 'gtest_main', 
                              ])
    objects_onlyfor = envTEST.Object(glob.glob('tests/*Tests.cpp'))
    program = envTEST.Program(testProg, objects_common + objects_onlyfor)

    envTESTDebug.Program(testProg + '.debug', debug_obj_common + makeObj(envTESTDebug, glob.glob('tests/*Tests.cpp'), '.debug.o'))
    # makeprogDebug(envDebug, glob.glob('tests/*Tests.cpp'), testProg+'.debug', debug_obj_common + )
    # program = makeprog(glob.glob('tests/*Tests.cpp'), testProg)
    # makeprogDebug(envTESTDebug, glob.glob('tests/*Tests.cpp'), testProg+'.debug.o', debug_obj_all)

    Command("tests.passed", testProg, runUnitTest)

    test_alias = Alias(testProg, [program], program[0].abspath)
    AlwaysBuild(test_alias)
else:
    print "Google test not found in /opt/local/lib/. Testing skipped."
####################################################
# def makeTags():
#     import subprocess
#     if os.path.exists('/opt/local/bin/gtags'):
#         subprocess.call(["/opt/local/bin/gtags", "-vv"])
# makeTags()
