#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/Utilities.h"
#include <list>
#include "cinder/Rand.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class particle {

public:
    particle();
    float angle, colorDist, currentColorDist, torque;
    Vec2f loc, dir, force;
    std::vector<Vec2f> interactLines;
    std::vector<float> interactStrengths;
    
};

particle::particle (){
    
}

class Particles_2DApp : public AppNative {

public:
    void prepareSettings(Settings* settings);
    void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
    void resize();
    
    // Drawing parameters
    float dotSize = 0, normalizedDotSize = 0.0065;
    float windowPadding = 20;
    Vec2f windowSize = Vec2f(600,800);
    Vec2f paddedWindowSize = windowSize+Vec2f(2*windowPadding,2*windowPadding);

    // Particle variables
    int numParticles = 700;
    int numPredators = 10;
    float preyVelocity = 0.003;
    std::vector<particle> particleVec, predatorVec;
    std::vector<particle>::iterator particleIter;
    std::vector<std::vector<Boolean>> inRangeFlags;
    std::vector<int> numberInRange, numberInColoringRange;
    std::vector<std::vector<float>> particleDistances;
    std::vector<float> interactOrientations;
    std::vector<Vec2f> interactionForces;
    
    // Interaction parameters
    float interactRange = 0, interactRangeNormalized = 0.035;
    float coloringRange = 10;
    float coloringScaling = 20;
    int meshSampling = 25;
    
    std::vector<std::vector<Vec2f>> binLocs;
    
    std::vector<std::vector<std::vector<particle*>>> binParticlePointers;
    std::vector<std::vector<int>> binNumParticles;
    std::vector<std::vector<std::vector<particle*>>> largeBinParticlePointers;
    std::vector<std::vector<int>> largeBinNumParticles;
    
    std::vector<std::vector<std::vector<particle*>>> binPredatorPointers;
    std::vector<std::vector<int>> binNumPredators;
    std::vector<std::vector<std::vector<particle*>>> largeBinPredatorPointers;
    std::vector<std::vector<int>> largeBinNumPredators;
    
    
    const float pi = 3.14159;
    double time = 0, deltat = 0;
    
};

void Particles_2DApp::prepareSettings(Settings* settings){
    
    settings->setFrameRate(70.0f);
    settings->setWindowSize(windowSize);
    
    if (windowSize[0]>windowSize[1]){
        interactRange = windowSize[1]*interactRangeNormalized;
        dotSize = windowSize[1]*normalizedDotSize;
    }else{
        interactRange = windowSize[0]*interactRangeNormalized;
        dotSize = windowSize[0]*normalizedDotSize;
    }
    
    paddedWindowSize = windowSize+Vec2f(2*windowPadding,2*windowPadding);
    
}

void Particles_2DApp::setup()
{
    
    // initialize particles
    particleVec.resize(0);
    for ( int kk = 0; kk < numParticles; kk++) {
        particleVec.push_back(particle());
        particleVec[kk].loc = Vec2f(randFloat(),randFloat());
        particleVec[kk].angle = randFloat();
        particleVec[kk].dir = Vec2f(cos(particleVec[kk].angle*2*pi),sin(particleVec[kk].angle*2*pi));
        particleVec[kk].currentColorDist = 100.0f;
        particleVec[kk].interactLines.resize(0);
        particleVec[kk].interactStrengths.resize(0);
    }
    
    // initialize predators
    predatorVec.resize(0);
    for ( int kk = 0; kk < numPredators; kk++) {
        predatorVec.push_back(particle());
        predatorVec[kk].loc = Vec2f(randFloat(),randFloat());
        predatorVec[kk].angle = randFloat();
        predatorVec[kk].dir = Vec2f(cos(predatorVec[kk].angle*2*pi),
                                    sin(predatorVec[kk].angle*2*pi));
        predatorVec[kk].currentColorDist = 100.0f;
        predatorVec[kk].interactLines.resize(0);
    }
    
    console() << predatorVec.size() << '\n';
    
    // Initialize vectors to handle particles
    std::vector<Boolean> inputVecBool;
    inputVecBool.resize(numParticles, false);
    inRangeFlags.resize(numParticles,inputVecBool);
    numberInRange.resize(numParticles,0);
    std::vector<float> inputVecDist;
    inputVecDist.resize(numParticles, 0.0f);
    particleDistances.resize(numParticles,inputVecDist);
    interactOrientations.resize(numParticles,0.0);
    interactionForces.resize(numParticles, Vec2f(0,0));
    
    
    
    for (int kk = 0;kk<meshSampling;kk++){
        std::vector<Vec2f> inputLocsBins;
        binLocs.push_back(inputLocsBins);
        for ( int ll = 0; ll<meshSampling; ll++){
            binLocs[kk].push_back(Vec2f((kk+0.5)*windowSize[0]/meshSampling,
                                        (ll+0.5)*windowSize[1]/meshSampling));
        }
    }
    
    
}

void Particles_2DApp::mouseDown( MouseEvent event )
{
    
    if( event.isLeft() && event.isShiftDown() ) {
        console() << "Full screen toggle" << std::endl;
        
        if (isFullScreen()) {
            setFullScreen(false);
            showCursor();
        }else{
            setFullScreen(true);
            hideCursor();
        }
            
    }
    
}

void Particles_2DApp::resize()
{
    windowSize = getWindowSize();
    paddedWindowSize = windowSize+Vec2f(2*windowPadding,2*windowPadding);
    
}

void Particles_2DApp::update()
{
    
    
    // Periodically change the interaction range
    
    
    float scalingFactor = (0.9+0.7*sin(time*2*pi/60.0));
    if (paddedWindowSize[0]>paddedWindowSize[1]){
        interactRange = paddedWindowSize[1]*interactRangeNormalized*scalingFactor;
        dotSize = paddedWindowSize[1]*normalizedDotSize*(0.3+0.7*scalingFactor);
    }else{
        interactRange = paddedWindowSize[0]*interactRangeNormalized*scalingFactor;
        dotSize = paddedWindowSize[0]*normalizedDotSize*(0.3+0.7*scalingFactor);
    }

    
    
    
    // Set up vectors and pointers for mesh binning
    std::vector<particle*> particlePointerVec;
    particlePointerVec.resize(0);
    std::vector<particle*> largeBinParticlePointerVec;
    largeBinParticlePointerVec.resize(0);
    std::vector<std::vector<particle*>> particlePointerVecVec,
    largeBinParticlePointerVecVec;
    particlePointerVecVec.resize(meshSampling,particlePointerVec);
    binParticlePointers.assign(meshSampling, particlePointerVecVec);
    largeBinParticlePointers.assign(meshSampling, particlePointerVecVec);

    // Set up vectors and pointers for predator mesh binning
    std::vector<particle*> predatorPointerVec;
    predatorPointerVec.resize(0);
    std::vector<particle*> largeBinPredatorPointerVec;
    largeBinPredatorPointerVec.resize(0);
    std::vector<std::vector<particle*>> predatorPointerVecVec,
    largeBinPredatorPointerVecVec;
    predatorPointerVecVec.resize(meshSampling,predatorPointerVec);
    binPredatorPointers.assign(meshSampling, predatorPointerVecVec);
    largeBinPredatorPointers.assign(meshSampling, predatorPointerVecVec);

    
    
    
    std::vector<int> inputNumVec;
    inputNumVec.resize(meshSampling,0);
    binNumParticles.assign(meshSampling,inputNumVec);
    largeBinNumParticles.assign(meshSampling,inputNumVec);
    
    binNumPredators.assign(meshSampling,inputNumVec);
    largeBinNumPredators.assign(meshSampling,inputNumVec);
    
    // --- Sort individual particles into mesh points

    int xBin, yBin;
    Vec2f thisLoc;
    
    for ( particleIter = particleVec.begin(); particleIter<particleVec.end(); particleIter++ ) {
        
        thisLoc = particleIter->loc;
        
        xBin = floor(thisLoc[0]*meshSampling);
        yBin = floor(thisLoc[1]*meshSampling);
        
        binParticlePointers[xBin][yBin].push_back(&(*particleIter));
        binNumParticles[xBin][yBin]++;
        particleIter->colorDist = 0.0f;
        
    }

    
    for ( particleIter = predatorVec.begin(); particleIter<predatorVec.end(); particleIter++ ) {
        
        thisLoc = particleIter->loc;
        
        xBin = floor(thisLoc[0]*meshSampling);
        yBin = floor(thisLoc[1]*meshSampling);
        
        binPredatorPointers[xBin][yBin].push_back(&(*particleIter));
        binNumPredators[xBin][yBin]++;
        particleIter->colorDist = 0.0f;
        
    }

    
    
    // --- Update particles as they were binned together ---
    
    int partNum = 0, largeBinPartNum = 0, predatorNum = 0, largeBinPredatorNum = 0;
    float partDist, angleDiff;
    Vec2f diffVec = Vec2f(0,0);
    
    for ( int mm = 0; mm<meshSampling; mm++){
        for ( int nn = 0; nn<meshSampling; nn++) {
            if (binNumParticles[mm][nn]>0){
                
                // --- work with the particles referenced by pointers
                particlePointerVec = binParticlePointers[mm][nn];
                partNum = binNumParticles[mm][nn];
                
                predatorPointerVec = binPredatorPointers[mm][nn];
                predatorNum = binNumPredators[mm][nn];
                
                // --- pool further away bins
                largeBinParticlePointerVec = binParticlePointers[mm][nn];
                largeBinPartNum = binNumParticles[mm][nn];

                largeBinPredatorPointerVec = binPredatorPointers[mm][nn];
                largeBinPredatorNum = binNumPredators[mm][nn];
                
                // add the neighboring boxes
                if (mm>0){
                    largeBinParticlePointerVec.insert(
                        largeBinParticlePointerVec.end(),
                        binParticlePointers[mm-1][nn].begin(),
                        binParticlePointers[mm-1][nn].end() );
                    largeBinPartNum += binNumParticles[mm-1][nn];
                    
                    largeBinPredatorPointerVec.insert(
                          largeBinPredatorPointerVec.end(),
                          binPredatorPointers[mm-1][nn].begin(),
                          binPredatorPointers[mm-1][nn].end() );
                    largeBinPredatorNum += binNumPredators[mm-1][nn];
                } else {
                    largeBinParticlePointerVec.insert(
                          largeBinParticlePointerVec.end(),
                          binParticlePointers[meshSampling-1][nn].begin(),
                          binParticlePointers[meshSampling-1][nn].end() );
                    largeBinPartNum += binNumParticles[meshSampling-1][nn];
                    
                    largeBinPredatorPointerVec.insert(
                          largeBinPredatorPointerVec.end(),
                          binPredatorPointers[meshSampling-1][nn].begin(),
                          binPredatorPointers[meshSampling-1][nn].end() );
                    largeBinPredatorNum += binNumPredators[meshSampling-1][nn];
                }
                
                if (nn>0){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm][nn-1].begin(),
                      binParticlePointers[mm][nn-1].end() );
                    largeBinPartNum += binNumParticles[mm][nn-1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm][nn-1].begin(),
                      binPredatorPointers[mm][nn-1].end() );
                    largeBinPredatorNum += binNumPredators[mm][nn-1];
                } else {
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm][meshSampling-1].begin(),
                      binParticlePointers[mm][meshSampling-1].end() );
                    largeBinPartNum += binNumParticles[mm][meshSampling-1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm][meshSampling-1].begin(),
                      binPredatorPointers[mm][meshSampling-1].end() );
                    largeBinPredatorNum += binNumPredators[mm][meshSampling-1];
                }

                if (mm<(meshSampling-1)){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm+1][nn].begin(),
                      binParticlePointers[mm+1][nn].end() );
                    largeBinPartNum += binNumParticles[mm+1][nn];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm+1][nn].begin(),
                      binPredatorPointers[mm+1][nn].end() );
                    largeBinPredatorNum += binNumPredators[mm+1][nn];
                } else {
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[0][nn].begin(),
                      binParticlePointers[0][nn].end() );
                    largeBinPartNum += binNumParticles[0][nn];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[0][nn].begin(),
                      binPredatorPointers[0][nn].end() );
                    largeBinPredatorNum += binNumPredators[0][nn];
                }

                if (nn<(meshSampling-1)){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm][nn+1].begin(),
                      binParticlePointers[mm][nn+1].end() );
                    largeBinPartNum += binNumParticles[mm][nn+1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm][nn+1].begin(),
                      binPredatorPointers[mm][nn+1].end() );
                    largeBinPredatorNum += binNumPredators[mm][nn+1];
                } else {
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm][0].begin(),
                      binParticlePointers[mm][0].end() );
                    largeBinPartNum += binNumParticles[mm][0];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm][0].begin(),
                      binPredatorPointers[mm][0].end() );
                    largeBinPredatorNum += binNumPredators[mm][0];
                }

                if (mm>0 & nn>0){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm-1][nn-1].begin(),
                      binParticlePointers[mm-1][nn-1].end() );
                    largeBinPartNum += binNumParticles[mm-1][nn-1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm-1][nn-1].begin(),
                      binPredatorPointers[mm-1][nn-1].end() );
                    largeBinPredatorNum += binNumPredators[mm-1][nn-1];
                }

                if (mm<(meshSampling-1) & nn>0){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm+1][nn-1].begin(),
                      binParticlePointers[mm+1][nn-1].end() );
                    largeBinPartNum += binNumParticles[mm+1][nn-1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm+1][nn-1].begin(),
                      binPredatorPointers[mm+1][nn-1].end() );
                    largeBinPredatorNum += binNumPredators[mm+1][nn-1];
                }

                if (mm>0 & nn<(meshSampling-1)){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm-1][nn+1].begin(),
                      binParticlePointers[mm-1][nn+1].end() );
                    largeBinPartNum += binNumParticles[mm-1][nn+1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm-1][nn+1].begin(),
                      binPredatorPointers[mm-1][nn+1].end() );
                    largeBinPredatorNum += binNumPredators[mm-1][nn+1];
                }
                
                if (mm<(meshSampling-1) & nn<(meshSampling-1)){
                    largeBinParticlePointerVec.insert(
                      largeBinParticlePointerVec.end(),
                      binParticlePointers[mm+1][nn+1].begin(),
                      binParticlePointers[mm+1][nn+1].end() );
                    largeBinPartNum += binNumParticles[mm+1][nn+1];
                    
                    largeBinPredatorPointerVec.insert(
                      largeBinPredatorPointerVec.end(),
                      binPredatorPointers[mm+1][nn+1].begin(),
                      binPredatorPointers[mm+1][nn+1].end() );
                    largeBinPredatorNum += binNumPredators[mm+1][nn+1];
                }
                
                
                numberInRange.assign(partNum,0);
                numberInColoringRange.assign(partNum,0);
                interactOrientations.assign(partNum,0.0);
                interactionForces.assign(partNum, Vec2f(0,0));
                

                // Get pairwise distances and interactions
                for ( int kk = 0; kk<partNum; kk++) {
                    
                    particlePointerVec[kk]->interactLines.resize(0);
                    particlePointerVec[kk]->interactStrengths.resize(0);
                    
                    for ( int ll = 0; ll<largeBinPartNum; ll++){

                        if (particlePointerVec[kk]!=largeBinParticlePointerVec[ll]) {

                            diffVec = (particlePointerVec[kk]->loc
                                       -largeBinParticlePointerVec[ll]->loc);
                            
                            if (diffVec[0]<-0.5) {
                                diffVec[0] += 1.0f;
                            }
                            if (diffVec[0]>0.5) {
                                diffVec[0] -= 1.0f;
                            }
                            if (diffVec[1]<-0.5) {
                                diffVec[1] += 1.0f;
                            }
                            if (diffVec[1]>0.5) {
                                diffVec[1] -= 1.0f;
                            }
                            diffVec *= paddedWindowSize;
                            partDist = pow(pow(diffVec[0],2.0)+pow(diffVec[1],2.0),0.5);
                            
                            if (partDist<=coloringRange){
                                
                                numberInColoringRange[kk] += 1;
                                particlePointerVec[kk]->colorDist += partDist;
                                
                            }
                            
                            if (partDist<=interactRange){
                                
                                numberInRange[kk] += 1;
                                
                                angleDiff = particlePointerVec[kk]->angle-largeBinParticlePointerVec[ll]->angle;
                                while (angleDiff<-0.5) {
                                    angleDiff += 1.0f;
                                }
                                while (angleDiff>0.5) {
                                    angleDiff -= 1.0f;
                                }

                                interactOrientations[kk] += angleDiff*(interactRange-partDist)/interactRange;
                                
                                interactionForces[kk] += diffVec*(interactRange-partDist)/interactRange
                                    +0.5*diffVec/(partDist+0.1);
                                
                            }
                            
                            if (partDist<=1.5*interactRange){
                                particlePointerVec[kk]->interactLines.push_back(diffVec);
                                particlePointerVec[kk]->interactStrengths.push_back(1.0-partDist/(1.5*interactRange));
                            }
                                
                        }
                        
                    }
                    
                    
                    // get predator interactions
                    
                    for ( int ll = 0; ll<largeBinPredatorNum; ll++) {
                        
                        diffVec = (particlePointerVec[kk]->loc
                                   -largeBinPredatorPointerVec[ll]->loc);
                        
                        if (diffVec[0]<-0.5) {
                            diffVec[0] += 1.0f;
                        }
                        if (diffVec[0]>0.5) {
                            diffVec[0] -= 1.0f;
                        }
                        if (diffVec[1]<-0.5) {
                            diffVec[1] += 1.0f;
                        }
                        if (diffVec[1]>0.5) {
                            diffVec[1] -= 1.0f;
                        }
                        diffVec *= paddedWindowSize;
                        partDist = pow(pow(diffVec[0],2.0)+pow(diffVec[1],2.0),0.5);
                        
                        if (partDist<=interactRange*2){
                            
                            interactionForces[kk] += 1*diffVec*(interactRange*2-partDist)/interactRange*2;
                            largeBinPredatorPointerVec[ll]->force += diffVec*(interactRange*2-partDist)/interactRange*2;
                            
                        }
                        
                    }
                    
                }
                
                
                // --- Update particle torque and force
                
                for (int kk = 0; kk<partNum; kk++) {
                    
                    // torque and coloring update
                    if (numberInRange[kk]>0){
                        particlePointerVec[kk]->torque = -interactOrientations[kk]/numberInRange[kk];
                        particlePointerVec[kk]->colorDist /= numberInColoringRange[kk];
                        particlePointerVec[kk]->currentColorDist +=
                            0.1f*(particlePointerVec[kk]->colorDist-particlePointerVec[kk]->currentColorDist);

                    } else {
                        particlePointerVec[kk]->torque = 0;
                        particlePointerVec[kk]->currentColorDist +=
                            0.01f*(100.0f-particlePointerVec[kk]->currentColorDist);
                    }
                    
                    // force update
                    particlePointerVec[kk]->force =
                            +preyVelocity*particlePointerVec[kk]->dir
                            +(0.1*interactionForces[kk])/paddedWindowSize;
                    

                }
                
            }
        }
    }

    
}

void Particles_2DApp::draw()
{

    // clear out the window with white
    gl::clear( Color( 0, 0, 0 ) );
    
    
    // Update particle orientation and position based on torque and force
    
    float newTime = ci::app::getElapsedSeconds();
    deltat = (newTime-time)*20;
    time = newTime;
    
    for ( int kk = 0; kk < numParticles; kk++) {
        
        // orientiation update
        particleVec[kk].angle += particleVec[kk].torque*deltat+0.05*(randFloat()-0.5)*sqrt(deltat);
        if (particleVec[kk].angle<0.0f) {
            particleVec[kk].angle += 1.0f;
        }
        if (particleVec[kk].angle>=1.0f) {
            particleVec[kk].angle -= 1.0f;
        }
        
        particleVec[kk].dir = Vec2f(cos(particleVec[kk].angle*2*pi),
                                    sin(particleVec[kk].angle*2*pi));
        
        // location update
        particleVec[kk].loc += particleVec[kk].force*deltat;
        
        while (particleVec[kk].loc[0]<0.0f) {
            particleVec[kk].loc[0] += 1.0f;
        }
        while (particleVec[kk].loc[0]>=1.0f) {
            particleVec[kk].loc[0] -= 1.0f;
        }
        
        while (particleVec[kk].loc[1]<0.0f) {
            particleVec[kk].loc[1] += 1.0f;
        }
        while (particleVec[kk].loc[1]>=1.0f) {
            particleVec[kk].loc[1] -= 1.0f;
        }
        
        
    }

    
    for ( int kk = 0; kk < numPredators; kk++) {
        
        // location update
        predatorVec[kk].loc += 0.15*predatorVec[kk].dir*deltat;
        predatorVec[kk].dir += 0.05*predatorVec[kk].force/paddedWindowSize*deltat - 0.2*predatorVec[kk].dir*deltat;
        predatorVec[kk].force -= predatorVec[kk].force*deltat;
        
        while (predatorVec[kk].loc[0]<0.0f) {
            predatorVec[kk].loc[0] += 1.0f;
        }
        while (predatorVec[kk].loc[0]>=1.0f) {
            predatorVec[kk].loc[0] -= 1.0f;
        }
        
        while (predatorVec[kk].loc[1]<0.0f) {
            predatorVec[kk].loc[1] += 1.0f;
        }
        while (predatorVec[kk].loc[1]>=1.0f) {
            predatorVec[kk].loc[1] -= 1.0f;
        }
        
        
    }

    
    
    
    // Do actual drawing

    Vec2f drawVector = Vec2f(0,0);

    gl::enableAlphaBlending();
    
//    for ( int kk = 0; kk < numPredators; kk++) {
//        
//        drawVector = predatorVec[kk].loc*paddedWindowSize;
//        
//        glColor4f(1,1,1,1.0);
//        
////        gl :: begin ( GL_POLYGON );
////        gl :: vertex ( drawVector[0]-dotSize*2 , drawVector[1]-dotSize*2 );
////        gl :: vertex ( drawVector[0]-dotSize*2 , drawVector[1]+dotSize*2 );
////        gl :: vertex ( drawVector[0]+dotSize*2 , drawVector[1]+dotSize*2 );
////        gl :: vertex ( drawVector[0]+dotSize*2 , drawVector[1]-dotSize*2 );
////        gl :: vertex ( drawVector[0]-dotSize*2 , drawVector[1]-dotSize*2 );
////        gl :: end ();
////        
//
//        gl::drawLine(drawVector-Vec2f(0,2*dotSize),drawVector+Vec2f(0,2*dotSize));
//        gl::drawLine(drawVector-Vec2f(2*dotSize,0),drawVector+Vec2f(2*dotSize,0));
//
//
//    }
    
    
    for ( int kk = 0; kk < numParticles; kk++) {
        
        drawVector = particleVec[kk].loc*paddedWindowSize-Vec2f(windowPadding,windowPadding);
        
        glColor4f(1,1,1,0.5);
        gl :: begin ( GL_POLYGON );
        gl :: vertex ( drawVector[0]-dotSize , drawVector[1]-dotSize );
        gl :: vertex ( drawVector[0]-dotSize , drawVector[1]+dotSize );
        gl :: vertex ( drawVector[0]+dotSize , drawVector[1]+dotSize );
        gl :: vertex ( drawVector[0]+dotSize , drawVector[1]-dotSize );
        gl :: vertex ( drawVector[0]-dotSize , drawVector[1]-dotSize );
        gl :: end ();
        
        int numInteract =  particleVec[kk].interactLines.size();
        
        if (numInteract>20){
            numInteract = 20;
        }
        
        glLineWidth(2.5f);
        if ( numInteract>0){
            for ( int nn = 0; nn<numInteract; nn++) {
                glColor4f(1,1,1,
                          particleVec[kk].interactStrengths[nn]);
                gl::drawLine(drawVector, drawVector-particleVec[kk].interactLines[nn]);
            }
        }
        
    }
    
}

CINDER_APP_NATIVE( Particles_2DApp, RendererGl )
