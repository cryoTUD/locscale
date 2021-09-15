# Run codes for proshade
import numpy as np
import proshade

def proshade_overlay(map1, map2, fitresol=4.0):
    pSet = proshade.ProSHADE_settings()
    pSet.task = proshade.OverlayMap
    pSet.verbose = -1
    pSet.setResolution(fitresol)
    pSet.firstModelOnly = False
    pSet.addStructure(map1)
    pSet.addStructure(map2)

    pRun = proshade.ProSHADE_run(pSet)  # Will take some time, maybe 30 seconds or so

    rotMatrix = proshade.getRotationMat(pRun)
    toOriginTranslation = proshade.getNumpyTranslationToOrigin(pRun)
    toMapCenTranslation = proshade.getNumpyTranslationToMapCentre(pRun)
    origToOverTranslation = proshade.getNumpyOriginToOverlayTranslation(pRun)

    rotmat = np.array(
        [
            [rotMatrix[0][0], rotMatrix[0][1], rotMatrix[0][2]],
            [rotMatrix[1][0], rotMatrix[1][1], rotMatrix[1][2]],
            [rotMatrix[2][0], rotMatrix[2][1], rotMatrix[2][2]],
        ],
        dtype="float",
    )
    return rotmat

def get_symmops_from_proshade(mapname):
    import sys
    import numpy
    import proshade

    ### Create the settings object
    pSet = proshade.ProSHADE_settings()
    ### Set settings values
    pSet.task = proshade.Symmetry
    pSet.verbose = 1
    pSet.setResolution(7.0)
    pSet.moveToCOM = False
    pSet.changeMapResolution = True
    pSet.changeMapResolutionTriLinear = False

    ### Create the structure object
    pStruct = proshade.ProSHADE_data(pSet)

    ### Read in the structure
    pStruct.readInStructure(mapname, 0, pSet)

    ### Process map
    pStruct.processInternalMap(pSet)

    ### Map to spheres
    pStruct.mapToSpheres(pSet)

    ### Compute spherical harmonics
    pStruct.computeSphericalHarmonics(pSet)

    ### Compute self-rotation function
    pStruct.computeRotationFunction(pSet)

    ### Detect symmetry
    #pStruct.detectSymmetryInStructurePython(pSet)
    pStruct.detectSymmetryInStructure(pSet)
    recSymmetryType = pStruct.getRecommendedSymmetryType(pSet)
    recSymmetryFold = pStruct.getRecommendedSymmetryFold(pSet)
    #recSymmetryAxes = proshade.getRecommendedSymmetryAxesPython(pStruct, pSet)
    #recSymmetryAxes = pStruct.getRecommendedSymmetryAxes(pSet) # works
    #recSymmetryAxes = pStruct.getAllCSyms(pSet)
    recSymmetryAxes = pStruct.getAllCSyms(pSet)
    #recSymmetryAxes = proshade.getSymmetryAxis(pStruct, pSet)

    ### Print results
    print("Detected " + str(recSymmetryType) + "-" + str(recSymmetryFold) + " symetry.")
    proshade_pg = str(recSymmetryType) + str(recSymmetryFold)
    print("Proshade point group: ", proshade_pg)
    print("Fold      x         y         z       Angle     Height")
    print(type(recSymmetryAxes))
    for iter in range(0, len(recSymmetryAxes)):
        print(
            "  %s    %+1.3f    %+1.3f    %+1.3f    %+1.3f    %+1.4f"
            % (
                recSymmetryAxes[iter][0],
                recSymmetryAxes[iter][1],
                recSymmetryAxes[iter][2],
                recSymmetryAxes[iter][3],
                recSymmetryAxes[iter][4],
                recSymmetryAxes[iter][5],
            )
        )
    fold, x, y, z, theta, peakh = [], [], [], [], [], []
    for row in recSymmetryAxes:
        fold.append(int(row[0]))
        x.append(row[1])
        y.append(row[2])
        z.append(row[3])
        theta.append(row[4])
        peakh.append(row[5])
    return [fold, x, y, z, peakh, proshade_pg]

