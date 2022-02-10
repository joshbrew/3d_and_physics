export class Math3D { //Bunch of static methods for doing math in 3D
    constructor() {

    }

	//Throwing a bunch in here for the hell of it
	static TWO_PI = Math.PI*2; //2PI
	static C = 299792458; //speed of light m/s
	static G = 6.67430e-11; //Newton's gravitation constant N*m^2 / kg^2
	static h = 6.62607015e-34; //Planck constant J*s
	static R = 8.31432e3; //Universal gas constant J / kg*mol*K
	static Ra = 287; //Air gas constant J / kg*K
	static H = 69.3; //Hubble constant km/s/Mpc 
	static kbar = 1.054571817e-34; //Dirac constant J*s
	static kB = 1.380649e-23; //Boltzmann constant J/K
	static ke = 8.9875517923e9; //Coulomb constant kg * m^3 * s^-2 * C^-2
	static me = 9.1093837015e-31; //electron mass kg
	static mp = 1.67262192369e-27; //proton mass kg
	static mn =	1.67492749804e-27; //neutron mass kg
	static P0 = 1.01325e5; //Sea level pressure N/m^2
	static T0 = 288.15; //Sea level room temperature K
	static p0 = 1.225; //Sea level air density kg/m^3
	static Na = 6.0220978e23; //Avogadro's number 1 / kg*mol
	static y = 1.405; //Adiabatic constant
	static M0 = 28.96643; //Sea level molecular weight
	static g0 = 9.80665; //Sea level gravity m/s^2
	static Re = 6.3781e6; //Earth radius m
	static B = 1.458e-6; //Thermal constant Kg / m*s*sqrt(kg)
	static S = 110.4; //Sutherland's constant K
	static Sigma = 3.65e-10; //Collision diameter of air m
	
    //2D integral approximation using rectangular area under the curve. If you need absolute values be sure to return that.
    static integral = (func=(x)=>{ let y=x; return y;}, range=[], stepx=0.01) => {
        let area = 0;
        for(let i = range[0]; i<range[1]; i+=stepx) {
            let y=func(i);
            area += y*stepx;
        }
        return area;
    }

    //3D double integral approximation
    static dintegral = (func=(x,y)=>{ let z = x+y; return z;}, range=[[],[]], stepx=0.01,stepy=stepx) => {
        let volume = 0;
        for(let i = range[0][0]+stepx; i<range[0][1]; i+=stepx) {
            for(let j = range[1][0]+stepy; j<range[1][1]; j+=stepy) {
                let z=func(i,j);
                volume += z*stepx*stepy;
            }
        }
        return volume;
    }

    //4D triple integral approximation
    static tintegral = (func=(x,y,z)=>{ let w=x+y+z; return w;}, range=[[],[],[]], stepx=0.01, stepy=stepx, stepz=stepx) => {
        let volume = 0;
        for(let i = range[0][0]+stepx; i<range[0][1]; i+=stepx) {
            for(let j = range[1][0]+stepy; j<range[1][1]; j+=stepy) {
                for(let k = range[2][0]+stepz; k<range[2][1]; k+=stepz) {
                    let w=func(i,j,k);
                    volume += w*stepx*stepy*stepz;
                }
            }
        }
        return volume;
    }

    //2D path integral approximation (the length of a curve)
    static pintegral = (func=(x)=>{ let y=x; return y; }, range=[], stepx=0.01) => {
        let length = 0;
        let y0 = undefined;
        let yi = undefined;
        for(let i = range[0]; i<range[1]; i+=stepx) {
            y0 = yi;
            yi = func(i);
            if(y0)
                length += this.distance([0,y0],[stepx,yi]);
        }
        return length;
    }

    static dot(vec1,vec2) { //nDimensional vector dot product
        let dot=0;
        for(let i=0; i<vec1.length; i++) {
            dot += vec1[i]*vec2[i];
        }
        return dot;
    }

    static cross3D(vec1,vec2) { //3D vector cross product
        return [
            vec1[1]*vec2[2]-vec1[2]*vec2[1],
            vec1[2]*vec2[0]-vec1[0]*vec2[2],
            vec1[0]*vec2[1]-vec1[1]*vec2[0]
        ];
    }

    static magnitude(vec) { //nDimensional magnitude
        let sqrd = 0;
        vec.forEach((c) => {
            sqrd+=c*c;
        })
        return Math.sqrt(sqrd)
    }

    static distance(point1, point2) { //nDimensional vector distance function
        let dsqrd = 0;
        point1.forEach((c,i) => {
            dsqrd += (point2[i] - c)*(point2[i] - c);
        })
        return Math.sqrt(dsqrd);
    }

    static makeVec(point1,point2) {  //Make vector from two nDimensional points (arrays)
        var vec = new Array(point1.length);
        point1.forEach((c,i) => {
            vec.push(point2[i]-c);
        })
        return vec;
    }

    static normalize(vec) { //nDimensional normalization
        let _mag = 1/this.magnitude(vec);
        let vecn = new Array(vec.length);
        vec.forEach((c,i) => {
            vecn[i] = c*_mag;
        })
        return vecn;
    }
    
    //arbitrary-sized vector operations
    static vecsub(vec,subvec) {
        let res = new Array(vec.length);
        for(let i = 0; i < vec.length; i++) {
            res[i] = vec[i] - subvec[i];
        }
        return res;
    } 

    static vecadd(vec1,vec2) {
        let res = new Array(vec1.length);
        for(let i = 0; i < vec1.length; i++) {
            res[i] = vec1[i] + vec2[i];
        }
        return res;
    } 

    static vecmul(vec1,vec2) {
        let res = new Array(vec1.length);
        for(let i = 0; i < vec1.length; i++) {
            res[i] = vec1[i] * vec2[i];
        }
        return res;
    } 
    
    static vecdiv(numvec,denvec) {
        let res = new Array(numvec.length);
        for(let i = 0; i < numvec.length; i++) {
            res[i] = numvec[i] / denvec[i];
        }
        return res;
    } 

    static vecscale(vec1,scalar) {
        let res = new Array(vec1.length);
        for(let i = 0; i < vec1.length; i++) {
            res[i] = vec1[i] * scalar;
        }
        return res;
    } 

    //Find normal to a plane define by points (v(1to2) cross v(1to3)), can set to return the reverse normal (v(1to3) cross v(1to2)). Use to calculate triangle normals
    static calcNormal(point1,point2,point3,pos=true) {
        var QR = makeVec(point1,point2);
        var QS = makeVec(point1,point3);

        if(pos === true){
            return Math3D.normalize(this.cross3D(QR,QS));
        }
        else {
            return Math3D.normalize(this.cross3D(QS,QR));
        }
    }

    //generate normals from a triangle poly mesh (assuming each triangle supplied follows the right hand rule)
    static calcNormalMesh(mesh=[1,2,3,4,5,6,7,8,9]){
        var normalMesh = new Array(mesh.length);
        for(var i = 0; i < mesh.length; i+=9) {
            var normal = this.calcNormal([mesh[i],mesh[i+1],mesh[i+2]],[mesh[i+3],mesh[i+4],mesh[i+5]],[mesh[i+6],mesh[i+7],mesh[i+8]]);
            normalMesh[ i ] = normal[0];
            normalMesh[i+1] = normal[1];
            normalMesh[i+2] = normal[2];
            normalMesh[i+3] = normal[0];
            normalMesh[i+4] = normal[1];
            normalMesh[i+5] = normal[2];
            normalMesh[i+6] = normal[0];
            normalMesh[i+7] = normal[1];
            normalMesh[i+8] = normal[2];
        }

        return normalMesh;
    }

    //Rotates a list of 3D vectors about the origin. Usually better to supply transforms as matrices for the GPU to multiply
    static rotateMesh(mesh=[1,2,3,4,5,6,7,8,9], pitch, roll, yaw) {
        var cosa = Math.cos(yaw);
        var sina = Math.sin(yaw);

        var cosb = Math.cos(pitch);
        var sinb = Math.sin(pitch);

        var cosc = Math.cos(roll);
        var sinc = Math.sin(roll);

        var Axx = cosa*cosb;
        var Axy = cosa*sinb*sinc - sina*cosc;
        var Axz = cosa*sinb*cosc + sina*sinc;

        var Ayx = sina*cosb;
        var Ayy = sina*sinb*sinc + cosa*cosc;
        var Ayz = sina*sinb*cosc - cosa*sinc;

        var Azx = -sinb;
        var Azy = cosb*sinc;
        var Azz = cosb*cosc;

        var result = [...mesh];

        for (var i = 0; i < mesh.length; i+=3) {
            var px = mesh[i];
            var py = mesh[i+1];
            var pz = mesh[i+2];

            result[i] = Axx*px + Axy*py + Axz*pz;
            result[i+1] = Ayx*px + Ayy*py + Ayz*pz;
            result[i+2] = Azx*px + Azy*py + Azz*pz;
        }

        return result;
    }

    //Mesh is an array of vec3's offset by idx+=3
    static translateMesh(mesh, xOffset, yOffset, zOffset) {
        var result = [...mesh];
        for(var i = 0; i < mesh.length; i+=3) {
            result[i]   = mesh[i]+xOffset
            result[i+1] = mesh[i+1]+yOffset;
            result[i+2] = mesh[i+2]+zOffset;
        }
    }

    //Scale about origin
    static scaleMesh(mesh, xScalar, yScalar, zScalar) {
        var result = [...mesh];
        for(var i = 0; i < mesh.length; i+=3) {
            result[i]   = mesh[i]*xScalar
            result[i+1] = mesh[i+1]*yScalar;
            result[i+2] = mesh[i+2]*zScalar;
        }
    }

    static transposeMat2D(mat2D){
		return mat2D[0].map((_, colIndex) => mat2.map(row => row[colIndex]));
    }

    static matmul2D(a, b) { //matmul2Dtiply two 2D matrices (array of arrays)
		var aNumRows = a.length, aNumCols = a[0].length,
			bNumRows = b.length, bNumCols = b[0].length,
			m = new Array(aNumRows);  // initialize array of rows
		for (var r = 0; r < aNumRows; ++r) {
		  m[r] = new Array(bNumCols); // initialize the current row
		  for (var c = 0; c < bNumCols; ++c) {
			m[r][c] = 0;             // initialize the current cell
			for (var i = 0; i < aNumCols; ++i) {
			  m[r][c] += a[r][i] * b[i][c];
			}
		  }
		}
		return m;
    }

    static makeIdentityM4() {
        return [
            [1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,1]
        ];
    }

    static makeTranslationM4(tx,ty,tz){
        return [
            [1,   0,  0, 0],
            [0,   1,  0, 0],
            [0,   0,  1, 0],
            [tx, ty, tz, 1]
        ];
    }

    static translateM4(mat4, tx, ty, tz) {
        var translate = this.makeTranslationM4(tx,ty,tz)

        return Math3D.matmul2D(mat4, translate);
    }

    static makeScaleM4(scaleX,scaleY,scaleZ){
        return [
            [scaleX, 0, 0, 0],
            [0, scaleY, 0, 0],
            [0, 0, scaleZ, 0],
            [0, 0,      0, 1]
        ];

    }

    static scaleM4(mat4,scaleX,scaleY,scaleZ){
        var scale = this.makeScaleM4(scaleX,scaleY,scaleZ);
        return Math3D.matmul2Dtiply(mat4, scale);
    }


    static xRotationM4(angleInRadians) {
        var c = Math.cos(angleInRadians);
        var s = Math.sin(angleInRadians);

        return [
          [1, 0, 0, 0],
          [0, c, s, 0],
          [0,-s, c, 0],
          [0, 0, 0, 1],
        ];
    }

    static yRotationM4(angleInRadians) {
        var c = Math.cos(angleInRadians);
        var s = Math.sin(angleInRadians);

        return [
          [c, 0,-s, 0],
          [0, 1, 0, 0],
          [s, 0, c, 0],
          [0, 0, 0, 1]
        ];
    }

    static zRotationM4(angleInRadians) {
        var c = Math.cos(angleInRadians);
        var s = Math.sin(angleInRadians);

        return [
           [c, s, 0, 0],
          [-s, c, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1]
        ];
    }

    static lookAtM4(source=[0,0,0], target=[1,1,1], up=[0,1,0]) {
        var zAxis = this.normalize([source[0]-target[0],source[1]-target[1],source[2]-target[2]]);
        var xAxis = this.normalize(this.cross3D(up, zAxis));
        var yAxis = this.normalize(this.cross3D(zAxis, xAxis));

        return [
            [ xAxis[0], xAxis[1], xAxis[2], 0],
            [ yAxis[0], yAxis[1], yAxis[2], 0],
            [ zAxis[0], zAxis[1], zAxis[2], 0],
            [source[0],source[1],source[2], 1]
        ];
    }

    //Rotate a 4D matrix
    static rotateM4(mat4, anglex, angley, anglez) {
        var result = [...mat4];
        if(anglex !== 0){
            result = Math3D.matmul2D(result,this.xRotationM4(anglex));
        }
        if(angley !== 0){
            result = Math3D.matmul2D(result,this.yRotationM4(angley));
        }
        if(anglez !== 0){
            result = Math3D.matmul2D(result,this.zRotationM4(anglez));
        }

        return result;
    }

    static rotatePoint1AboutPoint2(p1,p2,anglex,angley,anglez) {
        let rotatedM4 =
            Math3D.matmul2D(
                this.translateM4(
                    this.rotateM4(
                        this.makeTranslationM4(p2[0],p2[1],p2[2]),
                        anglex,angley,anglez),
                    -p2[0],-p2[1],-p2[2]),
                [...p1,1]
            );

        return [rotatedM4[0][3],rotatedM4[1][3],rotatedM4[2][3]]
    }

    //4D matrix inversion. This is atypical formatting (usually mat4s are represented by a 1D array, which is more efficient)
    static invertM4(mat4) {
        var m = mat4;
        var inv = [...mat4];
        inv[0][0] = m[1][1]  * m[2][2]* m[3][3] -
                m[1][1]  * m[2][3]* m[3][2]-
                m[2][1] * m[1][2] * m[3][3] +
                m[2][1] * m[1][3]* m[3][2]+
                m[3][1]* m[1][2] * m[2][3]-
                m[3][1]* m[1][3]* m[2][2];

        inv[1][0] = -m[1][0] * m[2][2]* m[3][3] +
                m[1][0] * m[2][3]* m[3][2]+
                m[2][0] * m[1][2] * m[3][3] -
                m[2][0] * m[1][3]* m[3][2]-
                m[3][0]* m[1][2] * m[2][3]+
                m[3][0]* m[1][3]* m[2][2];

        inv[2][0] = m[1][0] * m[2][1]* m[3][3] -
                m[1][0] * m[2][3]* m[3][1]-
                m[2][0] * m[1][1] * m[3][3] +
                m[2][0] * m[1][3]* m[3][1]+
                m[3][0]* m[1][1] * m[2][3]-
                m[3][0]* m[1][3]* m[2][1];

        inv[3][0] = -m[1][0] * m[2][1]* m[3][2]+
                m[1][0] * m[2][2]* m[3][1]+
                m[2][0] * m[1][1] * m[3][2]-
                m[2][0] * m[1][2]* m[3][1]-
                m[3][0]* m[1][1] * m[2][2]+
                m[3][0]* m[1][2]* m[2][1];

        inv[0][1] = -m[0][1] * m[2][2]* m[3][3] +
                m[0][1] * m[2][3]* m[3][2]+
                m[2][1] * m[0][2]* m[3][3] -
                m[2][1] * m[0][3]* m[3][2]-
                m[3][1]* m[0][2]* m[2][3]+
                m[3][1]* m[0][3]* m[2][2];

        inv[1][1] = m[0][0]  * m[2][2]* m[3][3] -
                m[0][0]  * m[2][3]* m[3][2]-
                m[2][0] * m[0][2]* m[3][3] +
                m[2][0] * m[0][3]* m[3][2]+
                m[3][0]* m[0][2]* m[2][3]-
                m[3][0]* m[0][3]* m[2][2];

        inv[2][1] = -m[0][0]  * m[2][1]* m[3][3] +
                m[0][0]  * m[2][3]* m[3][1]+
                m[2][0] * m[0][1]* m[3][3] -
                m[2][0] * m[0][3]* m[3][1]-
                m[3][0]* m[0][1]* m[2][3]+
                m[3][0]* m[0][3]* m[2][1];

        inv[3][1] = m[0][0]  * m[2][1]* m[3][2]-
                m[0][0]  * m[2][2]* m[3][1]-
                m[2][0] * m[0][1]* m[3][2]+
                m[2][0] * m[0][2]* m[3][1]+
                m[3][0]* m[0][1]* m[2][2]-
                m[3][0]* m[0][2]* m[2][1];

        inv[0][2] = m[0][1] * m[1][2]* m[3][3] -
                m[0][1] * m[1][3]* m[3][2]-
                m[1][1]  * m[0][2]* m[3][3] +
                m[1][1]  * m[0][3]* m[3][2]+
                m[3][1]* m[0][2]* m[1][3]-
                m[3][1]* m[0][3]* m[1][2];

        inv[1][2] = -m[0][0]  * m[1][2]* m[3][3] +
                m[0][0]  * m[1][3]* m[3][2]+
                m[1][0] * m[0][2]* m[3][3] -
                m[1][0] * m[0][3]* m[3][2]-
                m[3][0]* m[0][2]* m[1][3]+
                m[3][0]* m[0][3]* m[1][2];

        inv[2][2] = m[0][0]  * m[1][1] * m[3][3] -
                m[0][0]  * m[1][3]* m[3][1]-
                m[1][0] * m[0][1]* m[3][3] +
                m[1][0] * m[0][3]* m[3][1]+
                m[3][0]* m[0][1]* m[1][3]-
                m[3][0]* m[0][3]* m[1][1];

        inv[3][2] = -m[0][0]  * m[1][1] * m[3][2]+
                m[0][0]  * m[1][2]* m[3][1]+
                m[1][0] * m[0][1]* m[3][2]-
                m[1][0] * m[0][2]* m[3][1]-
                m[3][0]* m[0][1]* m[1][2]+
                m[3][0]* m[0][2]* m[1][1];

        inv[0][3] = -m[0][1]* m[1][2]* m[2][3]+
                m[0][1]* m[1][3]* m[2][2]+
                m[1][1] * m[0][2]* m[2][3]-
                m[1][1] * m[0][3]* m[2][2]-
                m[2][1]* m[0][2]* m[1][3]+
                m[2][1]* m[0][3]* m[1][2];

        inv[1][3] = m[0][0] * m[1][2]* m[2][3]-
                m[0][0] * m[1][3]* m[2][2]-
                m[1][0]* m[0][2]* m[2][3]+
                m[1][0]* m[0][3]* m[2][2]+
                m[2][0]* m[0][2]* m[1][3]-
                m[2][0]* m[0][3]* m[1][2];

        inv[2][3] = -m[0][0] * m[1][1] * m[2][3]+
                m[0][0] * m[1][3]* m[2][1]+
                m[1][0]* m[0][1]* m[2][3]-
                m[1][0]* m[0][3]* m[2][1]-
                m[2][0]* m[0][1]* m[1][3]+
                m[2][0]* m[0][3]* m[1][1];

        inv[3][3] = m[0][0] * m[1][1] * m[2][2]-
                m[0][0] * m[1][2]* m[2][1]-
                m[1][0]* m[0][1]* m[2][2]+
                m[1][0]* m[0][2]* m[2][1]+
                m[2][0]* m[0][1]* m[1][2]-
                m[2][0]* m[0][2]* m[1][1];

        return inv;
    }

    //Convert 2D matrix (array of arrays) to a Float32Array buffer
    static bufferMat2D(mat2D){
        var arraybuffer = [];
        mat2D.forEach((arr,i)=>{
            arraybuffer.push(...arr);
        });

        return new Float32Array(arraybuffer);
    }

    //Fairly efficient nearest neighbor search. Supply list of coordinates (array of Array(3)) and maximum radius to be considered a neighbor.
    //Returns a list of nodes with [{idx:0,neighbors:[{idx:j,position:[x,y,z],dist:d}]},{...},...]. Neighbors are auto sorted by distance.
    //Current complexity: n(n+1)/2, there are faster ways to do it but this should be good enough
    static nearestNeighborSearch(positions, isWithinRadius) {

        let node = {
            idx: null,
            position: [0,0,0],
            neighbors: []
        };

        let neighbor = {
            idx: null,
            position: [0,0,0],
            dist: null
        };

        var tree = [];

        for(var i = 0; i < positions.length; i++){
            let newnode = JSON.parse(JSON.stringify(node));
            newnode.idx = i;
            newnode.position = positions[i];
            tree.push(newnode);
        }

        //Nearest neighbor search. This can be heavily optimized.
        for(var i = 0; i < tree.length; i++) { //for each node
            for(var j = i+1; j < tree.length; j++) { //for each position left to check
                var dist = Math3D.distance(tree[i].position,tree[j].position);
                if(dist < isWithinRadius){
                    var newNeighbori = JSON.parse(JSON.stringify(neighbor));
                    newNeighbori.position = positions[j];
                    newNeighbori.dist = dist;
                    newNeighbori.idx = tree[j].idx;
                    tree[i].neighbors.push(newNeighbori);
                    var newNeighborj = JSON.parse(JSON.stringify(neighbor)); //push corresponding neighbor
                    newNeighborj.position = positions[i];
                    newNeighborj.dist = dist;
                    newNeighborj.idx = tree[i].idx;
                    tree[j].neighbors.push(newNeighborj);
                }
            }
            tree[i].neighbors.sort(function(a,b) {return a.dist - b.dist}); //Sort by distance, nearest to farthest
        }

        return tree;
    }

    //fibonacci sphere mesh
    static FibSphere(nPoints){
        var goldenRatio = (1 + Math.sqrt(5)) * .5;
        var goldenAngle = (2.0 - goldenRatio) * (2.0*Math.PI);

        var vertices = [];

        for(var i = 0; i<nPoints; i++){
            var t = i/nPoints;
            var angle1 = Math.acos(1-2*t);
            var angle2 = goldenAngle*i;

            var x = Math.sin(angle1)*Math.cos(angle2);
            var y = Math.sin(angle1)*Math.sin(angle2);
            var z = Math.cos(angle1);

            vertices.push(x,y,z);
        }

        return vertices; // Returns vertex list [x0,y0,z0,x1,y1,z1,...]

    }

}