import {Math3D} from './Math3D'

export class Camera { //pinhole camera model. Use to set your 3D rendering view model (untested)
    constructor (
        position=[0,0,0],
        target=[0,100,0],
        up=[0,1,0],
        clientWidth=window.innerWidth,
        clientHeight=window.innerHeight
        )
        {

        this.position = {x:position[0],y:position[1],z:position[2]};
        this.target = {x:target[0],y:target[1],z:target[2]};
        this.up = {x:up[0],y:up[1],z:up[2]};

        this.fov = 90;
        this.aspect = clientWidth/clientHeight;

        //View distance
        this.near = 0;
        this.far = 1000;

        //Focal length
        this.fx = 1;
        this.fy = 1;

        //Center image pixel location?
        this.cx = clientWidth*.5;
        this.cy = clientHeight*.5;

        this.cameraMat = this.getLookAtViewProjectionMatrix(position,target,up);

    }

    getPerspectiveMatrix(fieldOfViewInRadians=this.fov, aspectRatio=this.aspect, near=this.near, far=this.far) {
        var f = 1.0 / Math.tan(fieldOfViewInRadians / 2);
        var rangeInv = 1 / (near - far);

        return [
          [f / aspectRatio, 0,                          0,   0],
          [0,               f,                          0,   0],
          [0,               0,    (near + far) * rangeInv,  -1],
          [0,               0,  near * far * rangeInv * 2,   0]
        ];
    }


    getProjectionMatrix(width, height, depth) {
        return [
            [2/width,   0,        0, 0],
            [0, -2/height,        0, 0],
            [0,         0,  2/depth, 0],
            [-1,        1,        0, 1]
        ];
    }


    getCameraMatrix(fx=this.fx, fy=this.fy, cx=this.cx, cy=this.cx) {
        return [
            [fx, 0, cx, 0],
            [0, fy, cy, 0],
            [0,  0,  1, 0]
        ];
    }

    getCameraViewProjectionMatrix(xPos=this.position.x,yPos=this.position.y,zPos=this.position.z,rotX=0,rotY=0,rotZ=0) { //Translate geometry based on this result then set location. Demonstrated: https://webglfundamentals.org/webgl/lessons/webgl-3d-camera.html
        var cameraMat = this.getCameraMatrix(this.fx,this.fy);
        cameraMat = Math3D.rotateM4(cameraMat,rotX,rotY,rotZ);
        cameraMat = Math3D.translateM4(cameraMat, xPos, yPos, zPos);
        var viewMat = Math3D.invertM4(cameraMat);
        return Math3D.matmul2D(this.getPerspectiveMatrix(), viewMat); //View projection matrix result
    }

    getLookAtViewProjectionMatrix(source, target, up) {
        var cameraMat = Math3D.lookAtM4(source,target,up);
        var viewMat = Math3D.invertM4(cameraMat);
        return Math3D.matmul2D(this.getPerspectiveMatrix(), viewMat);
    }

    getCameraTransform() {
        return Float32Array(Math3D.bufferMat2D(this.cameraMat));
    }

    updateRotation(xRot=0,yRot=0,zRot=0) {
        this.cameraMat = Math3D.rotateM4(this.cameraMat, xRot, yRot, zRot);
    }

    updateTranslation(xPos=this.position.x,yPos=this.position.y,zPos=this.position.z) {
        this.cameraMat = Math3D.translateM4(this.cameraMat, xPos, yPos, zPos);
    }

    //Rotation in Radians
    rotateCameraAboutPoint(xPos=0,yPos=0,zPos=0,xRot,yRot,zRot){
        var anchorPoint = [xPos,yPos,zPos];
        var rotatedPosition = Math3D.rotatePoint1AboutPoint2(this.position,anchorPoint,xRot,yRot,zRot);
        this.position = {x:rotatedPosition[0],y:rotatedPosition[1],z:rotatedPosition[2]};
        this.updateTranslation();
        this.updateRotation(-xRot,-yRot,-zRot);
    }

    moveCamera(xPos = 0, yPos = 0, zPos = 0) {
        this.position = {x:xPos, y:yPos, z:zPos};
        this.updateTranslation();
    }

}
