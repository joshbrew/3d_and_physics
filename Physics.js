import {Math3D} from './Math3D'

//untested physics body model
//contains basic physics bodies with sphere, box, and point collisions (and the combinations)
//dynamic AABB and Octrees
// gravitational bodies
// closest point on line/polygon or is point inside box
export class Physics {
    constructor(nBodies = 10) {

        this.physicsBodies = [];
        this.boundingTree = [];

        this.globalSettings = {
            maxDistCheck: 1000,
            gravity: 9.81
        };

        this.dynamicBoundingVolumeTree = {
            proto:{
                parent:undefined,
                children:[],
                bodies:[],
                collisionType:"Box",
                collisionRadius: 1, 
                collisionBoundsScale: [1,1,1], //radius of bounding box
                position:{x:0,y:0,z:0} //center of bounds
            },
            tree:{} //head will be a hard copy of the prototype and so on
        };

        this.bodySettings = {

            index: null,

            collisionEnabled: true,
            collisionType: "Sphere", //Sphere, Box, Point
            collisionRadius: 1, //Radius of sphere or nearest point on side planes in a box
            collisionBoundsScale: [1,1,1], //Can distort the bounding box dimensions, doesn't affect the sphere yet.

            dynamic: true, //Does this object move? Changes objects added to the dynamic bounding volume

            position: [0,0,0], //[x,y,z] or [i,j,k]
            rotation: [0,0,0],
            velocity: [0,0,0],
            acceleration: [0,0,0],
            forceImpulse: [0,0,0], //Instantaneous force (resets after applying)

            mass: 1,
            drag: 0,
            restitution: 1, //Bounciness
            friction: 0, //Amount this surface slows other objects in contact along the contact plane

            attractor: false, //N-body attractor
            attractionAccel: 9.81,

            trigger: false,
            triggerFunc: null,

            parent:null,
            children: null, //Child object class instance (for modifying parameters)
        }

        for (let i = 0; i < nBodies; i++) {
            this.physicsBodies.push(JSON.parse(JSON.stringify(this.bodySettings)));
            this.physicsBodies[i].index = i;
        }
    }

    addBody(settings={},bodies=this.physicsBodies) {
        let newbody = JSON.parse(JSON.stringify(this.bodySettings));
        Object.assign(newbody, settings); //overwrite e.g. position
        bodies[i].index = bodies.length;
        bodies.push(newbody); 

        return newbody;
    }

    //set any values on the object
    setBody(body,settings={}) {
        Object.assign(body,settings);
        return true;
    }

    removeBody(index,bodies=this.physicsBodies) {
        if(!index) return false;
        if(index > bodies.length - 1) return false;
        bodies.splice(index,1);
        bodies.map((o) => {if(o.index > index) o.index--;});

        return true;
    }

    newBoundingVolume(settings={}) {
        let newvolume = JSON.parse(JSON.stringify(this.dynamicBoundingVolumeTree.proto));
        Object.assign(newvolume, settings);

        return newvolume;
    }

    //dynamic AABB trees or octrees
    generateBoundingVolumeTree(bodies=this.physicsBodies, octree=true) {
        
        /*
        How to make dynamic bounding volume tree:
        1. Find the bounding volume of all of the objects combined.
        2. Now subdivide bounding volumes by whatever magnitude until you have groups of 2-5 objects close to each other.
        3. Use box collision checks on the tree to find the best candidates to search for collisions
        */

        let maxX, maxY, maxZ;
        let minX = 0, minY = 0, minZ = 0;
        let positions = [];
        let minRadius = this.globalSettings.maxDistCheck;

        bodies.forEach((body)=>{

            let xx = body.position[0]+body.collisionRadius*body.collisionBoundsScale[0];
            let yy = body.position[1]+body.collisionRadius*body.collisionBoundsScale[1];
            let zz = body.position[2]+body.collisionRadius*body.collisionBoundsScale[2];

            if(maxX < xx) maxX = xx;
            if(minX > xx) minX = xx;
            if(maxY < yy) maxY = yy;
            if(minY > yy) minY = yy;
            if(maxZ < zz) maxZ = zz;
            if(minZ > zz) minZ = zz;

            if(minRadius > body.collisionRadius) minRadius = body.collisionRadius;

            positions.push(body.position[0],body.position[1],body.position[2]);

        });

        let head = JSON.parse(JSON.stringify(this.dynamicBoundingVolumeTree.proto));

        let boxdims = [maxX-minX,maxY-minY,maxZ-minZ]; //height, width, depth 
        let boxbounds = [boxdims[0]*0.5,boxdims[1]*0.5,boxdims[2]*0.5]
        let boxpos = {x:boxbounds[0]+minX,y:boxbounds[1]+minY,z:boxbounds[2]+minZ}; //center of box
        
        head.position = boxpos;
        head.collisionBoundsScale = boxbounds; //radius to centers of each sides i.e. distance from center

        head.bodies = bodies;
        
        this.dynamicBoundingVolumeTree.tree = head;

        minRadius *= 2;

        if(octree === true) { //octrees
            function genOct(parentPos,halfbounds) { //return center positions of child octree cubes, radius = parent half radius
                let oct1 = [parentPos[0]+halfbounds[0],parentPos[1]+halfbounds[1],parentPos[2]+halfbounds[2]]; //+x+y+z
                let oct2 = [parentPos[0]-halfbounds[0],parentPos[1]+halfbounds[1],parentPos[2]+halfbounds[2]]; //-x+y+z
                let oct3 = [parentPos[0]+halfbounds[0],parentPos[1]-halfbounds[1],parentPos[2]+halfbounds[2]]; //+x-y+z
                let oct4 = [parentPos[0]+halfbounds[0],parentPos[1]+halfbounds[1],parentPos[2]-halfbounds[2]]; //+x+y-z
                let oct5 = [parentPos[0]-halfbounds[0],parentPos[1]-halfbounds[1],parentPos[2]+halfbounds[2]]; //-x-y+z
                let oct6 = [parentPos[0]-halfbounds[0],parentPos[1]+halfbounds[1],parentPos[2]-halfbounds[2]]; //-x+y-z
                let oct7 = [parentPos[0]+halfbounds[0],parentPos[1]-halfbounds[1],parentPos[2]-halfbounds[2]]; //+x-y-z
                let oct8 = [parentPos[0]-halfbounds[0],parentPos[1]-halfbounds[1],parentPos[2]-halfbounds[2]]; //-x-y-z

                return [oct1,oct2,oct3,oct4,oct5,oct6,oct7,oct8,oct9];
            }

            function genOctTree(head) {           
                let halfbounds = [head.collisionBoundsScale[0]*0.5,head.collisionBoundsScale[1]*0.5,head.collisionBoundsScale[2]*0.5];     
                let octPos = genOct(head.position,halfbounds);
                let check = [...head.bodies];
                for(let i = 0; i < 8; i++) {
                    let newvolume = this.newBoundingVolume({position:octPos[i],collisionBoundsScale:halfbounds});
                    newvolume.parent = head;
                    //now check if any of the bodies are within these and eliminate from the check array
                    for(let j = check.length-1; j >= 0; j++) {
                        let collided = this.collisionCheck(check[j],newvolume);
                        if(collided) {
                            newvolume.bodies.push(check[j]);
                            check.splice(j,1);
                        }
                    } //recursively check each oct for bodies until only 3 bodies are contained
                    if(newvolume.bodies > 2) {
                        head.children.push(newvolume);
                        newvolume.parent = head;
                        if(newvolume.bodies > 3 && newvolume.collisionRadius*0.5 > minRadius) {genOctTree(newvolume);}
                    }
                }
            }

            genOctTree(head);
            
            return head;
        }
        else { //dynamic AABB trees
            let tree = Math3D.nearestNeighborSearch(positions,this.globalSettings.maxDistCheck);

            /**
            let node = { //all positions
                idx: null, //index in input position array
                position: [0,0,0], //copy of position vector
                neighbors: [] //neighbor nodes within the limit distance
            };

            let neighbor = {
                idx: null, //index in input position array
                position: [0,0,0], //copy of position vector
                dist: null //distance from parent
            };
             */

            let index = Math.floor(Math.random()*tree.length); //beginning with random node
            let searching = true; 
            let count = 0;

            let genBoundingBoxLevel = (tree,volumes) => {
                let newVolumes = [];
                let volumePositions = [];
                let foundidxs=[];
                while(searching && count < tree.length) { 
                    let node = tree[index]; 
                    let i = 0; //starting with the farthest neighbor in the node
                    let j = 0;
                    let ux = 0, uy = 0, uz = 0, mx = 0, my = 0, mz = 0;
                    let newvolume = this.newBoundingVolume()
                    while(i < node.neighbors.length && j < 3) { //make a box around the nearest 3 neighbors 
                        if(foundidxs.indexOf(node.neighbors[i].idx) > -1) { i++; continue; }
                        if(ux > positions[node.neighbors[i].idx][0]-volumes[node.neighbors[i].idx].collisionBoundsScale[0]) ux = positions[node.neighbors[i].idx][0]-volumes[node.neighbors[i].idx].collisionBoundsScale[0];
                        if(mx < positions[node.neighbors[i].idx][0]+volumes[node.neighbors[i].idx].collisionBoundsScale[0]) mx = positions[node.neighbors[i].idx][0]+volumes[node.neighbors[i].idx].collisionBoundsScale[0];
                        if(uy > positions[node.neighbors[i].idx][1]-volumes[node.neighbors[i].idx].collisionBoundsScale[1]) uy = positions[node.neighbors[i].idx][1]-volumes[node.neighbors[i].idx].collisionBoundsScale[1];
                        if(my < positions[node.neighbors[i].idx][1]+volumes[node.neighbors[i].idx].collisionBoundsScale[1]) my = positions[node.neighbors[i].idx][1]+volumes[node.neighbors[i].idx].collisionBoundsScale[1];
                        if(uz > positions[node.neighbors[i].idx][2]-volumes[node.neighbors[i].idx].collisionBoundsScale[2]) uz = positions[node.neighbors[i].idx][2]-volumes[node.neighbors[i].idx].collisionBoundsScale[2];
                        if(mz < positions[node.neighbors[i].idx][2]+volumes[node.neighbors[i].idx].collisionBoundsScale[2]) mz = positions[node.neighbors[i].idx][2]+volumes[node.neighbors[i].idx].collisionBoundsScale[2];
                        newvolume.children.push(volumes[node.neighbors[i].idx]);
                        volumes[node.neighbors[i].idx].parent = newvolume;
                        tree.splice(node.neighbors[i].idx,1);
                        foundidxs.push(node.neighbors[i].idx);
                        i++; j++;
                    }
                    let bounds = [mx*0.5,my*0.5,mz*0.5];
                    let pos = [bounds[0]+ux,bounds[1]+uy,bounds[2]+uz]
                    newvolume.position = pos;
                    newvolume.collisionBoundsScale = bounds;
                    
                    //then walk to the nearest unchecked node and make a box around the next 2 or 3 nodes
                    // then continue till you run out of nodes to check. Shoud be left with a bounding tree with larger to smaller boxes
                    newVolumes.push(newvolume);
                    volumePositions.push(pos);
                    
                    //now find the next not-found neighbor
                    while(node.neighbors.length > i) {
                        if(foundidxs.indexOf(node.neighbors[i].idx) < 0) break;
                            i++;
                    }
                    if(node.neighbors.length > i) {
                        index = node.neighbors.length[i].idx;
                    } else index = 0; //just jump back to zero and keep looking
                    
                    count++;
                }

                return [newVolumes,volumePositions];
            }

            let result = genBoundingBoxLevel(tree,volumes);
            let nextTree = Math3D.nearestNeighborSearch(result[1],this.globalSettings.maxDistCheck);
            while(nextTree.length > 2) {
                result = genBoundingBoxLevel(nextTree,result[0]);
                nextTree = Math3D.nearestNeighborSearch(result[1],this.globalSettings.maxDistCheck);
            }

            head.children = result[0]; //that should parent the final bounding boxes to the main box

            return head;
        }
        /*
            Now search children recursively till you have small enough groups to speed up collision checking

            nearestneighborsearch bodies by radius,
            next create bounds at that radius to capture all bodies
            repeat at a smaller radius till no more objects can be grouped at a minimum group size or radius
            bodies will be referenced in each 
        */
    }

    timeStep(dt, bodies=this.physicsBodies) { //dt in seconds

        /* //Nearest neighbor search optimization for collision detection (to cut down array searching), can make this only fire every n-milliseconds for speed
        var neighborTree = Math3D.nearestNeighborSearch(positions,this.globalSettings.maxDistCheck);
        neighborTree.forEach((node,i) => {
            var body1 = this.physicsBodies[i];
            node.neighbors.forEach((neighbor,j) => {
                var body2 = this.physicsBodies[j];
                var isColliding = this.collisionCheck(body,otherBody);
                if(isColliding === true) { resolveCollision(body,otherBody); }
            });
        });
        */

        bodies.forEach((body,i) => {

            //var positions = new Array(this.physicsBodies.length);

            for(var j = i+1; j < bodies.length; j++) {
                var otherBody = bodies[j];

                //Collision Check
                var isColliding = this.collisionCheck(body,otherBody);
                if(isColliding === true) {
                    this.resolveCollision(body,otherBody); //Now calculate forces
                    this.resolveCollision(otherBody,body); //Now calculate forces
                }

                //Gravity check
                if(body.attractor === true && otherBody.attractor === true) {
                    this.resolveAttractor(body,otherBody);
                }

            }

            //Resolve Attractors

            //Apply any forces
            body.acceleration[0] += forceImpulse[0]/body.mass - body.acceleration[0]*drag;
            body.acceleration[0] += forceImpulse[1]/body.mass - body.acceleration[1]*drag;
            body.acceleration[0] += forceImpulse[2]/body.mass - body.acceleration[2]*drag - this.globalSettings.gravity*dt;

            body.forceImpulse[0] = 0;
            body.forceImpulse[1] = 0;
            body.forceImpulse[2] = 0;

            body.velocity[0] += body.acceleration[0]*dt;
            body.velocity[1] += body.acceleration[1]*dt;
            body.velocity[2] += body.acceleration[2]*dt;

            //Finally, calculate new positions
            body.position[0] += body.velocity[0]*dt;
            body.position[1] += body.velocity[1]*dt;
            body.position[2] += body.velocity[2]*dt;
        });
    }

    // V = Vold*dt + a*dt^2, basic projectile motion equation
    calcVelocityStep(vOld=[0,0,0],a=[0,0,0],dt) {
        return [
            vOld[0]*dt + a[0]*dt*dt,
            vOld[1]*dt + a[1]*dt*dt,
            vOld[2]*dt + a[2]*dt*dt
        ];
    }

    // F = m*a for 3D vecs
    calcForce(m, a=[0,0,0]) {
        return [
            m*a[0],
            m*a[1],
            m*a[2]
        ];
    }

    // a = F/m for 3D vecs
    calcAccelFromForce(F=[0,0,0], m=0) {
        return [
            F[0]/m,
            F[1]/m,
            F[2]/m
        ];
    }

    resolveCollision(body1,body2) { //Resolve what body1 does in contact with body2 (call twice with bodies reversed to calculate in both directions)
        //Elastic collisions
        var directionVec = Math3D.makeVec(body1.position,body2.position); //Get direction toward body2
        var normal = Math3D.normalize(directionVec);
        if(body2.collisionType === "Sphere" || body2.collisionType === "Point") {

            var body1velocityMag = Math3D.magnitude(body1.velocity);

            var body2AccelMag = Math3D.magnitude(body2.acceleration);
            var body2AccelNormal = Math3D.normalize(body2.acceleration);

            body1.velocity = [-normal[0]*body1velocityMag*body1.restitution,-normal[1]*body1velocityMag*body1.restitution,-normal[2]*body1velocityMag*body1.restitution]; //Adjust velocity

            body1.forceImpulse[0] -= body2AccelMag*body2AccelNormal[0]*body2.mass;
            body1.forceImpulse[1] -= body2AccelMag*body2AccelNormal[1]*body2.mass;
            body1.forceImpulse[2] -= body2AccelMag*body2AccelNormal[2]*body2.mass;

        }
        else if (body2.collisionType === "Box") {
            //Find which side was collided with
            var max = Math.max(...directionVec);
            var min = Math.min(...directionVec);
            var side = max;;
            if(Math.abs(min) > max) {
                side = min;
            }
            var idx = directionVec.indexOf(side);

            body1.velocity[idx] = -body1.velocity[idx]*body1.restitution; //Reverse velocity

            var body2AccelMag = Math3D.magnitude(body2.acceleration);
            var body2AccelNormal = Math3D.normalize(body2.acceleration);

            body1.forceImpulse[idx] = -body2AccelNormal[idx]*body2AccelMag*body2.mass; //Add force

            //Apply Friction
        }
        return;
    };

    resolveAttractor(body1,body2) {
        if(body1.mass == 0 && body2.mass === 0) return;
        //Gravitational pull of nBodies
        var dist = Math3D.distance(body1.position,body2.position);
        var vec1 = Math3D.normalize(Math3D.makeVec(body1.position,body2.position)); // a to b
        var vec2 = Math3D.normalize(Math3D.makeVec(body2.position,body1.position)); // b to a

        //Newton's law of gravitation
        var Fg = 0.00000000006674 * body1.mass * body2.mass / (dist*dist);

        //Get force vectors
        FgOnBody1 = [vec1[0]*Fg,vec1[1]*Fg,vec1[2]*Fg];
        FgOnBody2 = [vec2[0]*Fg,vec2[1]*Fg,vec2[2]*Fg];

        body1.forceImpulse[0] += FgOnBody1[0];
        body1.forceImpulse[1] += FgOnBody1[1];
        body1.forceImpulse[2] += FgOnBody1[2];

        body2.forceImpulse[0] += FgOnBody2[0];
        body2.forceImpulse[1] += FgOnBody2[1];
        body2.forceImpulse[2] += FgOnBody2[2];
        return;
    }

    //Checks if two bodies are colliding based on their collision setting
    collisionCheck(body1,body2) {
        if(body1.collisionEnabled === false || body2.collisionEnabled === false) return false;

        //Check if within a range close enough to merit a collision check
        if(Math3D.distance(body1.position,body2.position) < Math.max(...body1.scale)*body1.collisionRadius+Math.max(...body2.scale)*body2.collisionRadius) {
            //Do collision check
            let isColliding = false;
            if(body1.collisionType === "Sphere") {
                if(body2.collisionType === "Sphere") { isColliding = this.sphericalCollisionCheck(body1,body2);}
                if(body2.collisionType === "Box") { isColliding = this.sphereBoxCollisionCheck(body1,body2);}
                if(body2.collisionType === "Point") { isColliding = this.isPointInsideSphere(body2.position,body1);}
            }
            else if(body1.collisionType === "Box" ) {
                if(body2.collisionType === "Sphere") { isColliding = this.sphereBoxCollisionCheck(body2,body1);}
                if(body2.collisionType === "Box") { isColliding = this.boxCollisionCheck(body1,body2);}
                if(body2.collisionType === "Point") { isColliding = this.isPointInsideBox(body1.position,body1); }
            }
            else if (body1.collisionType === "Point") {
                if(body2.collisionType === "Sphere") { isColliding = this.isPointInsideSphere(body1.position,body2); }
                if(body2.collisionType === "Box") { isColliding = this.isPointInsideBox(body1.position,body2); }
            }

            return isColliding;
        }
        else return false


    }

    //Check if point is inside the spherical volume
    isPointInsideSphere(point=[0,0,0],body) {
        let dist = Math3D.distance(point,body.position);

        return dist < body.collisionRadius;
    }

    //Collision between two spheres
    sphericalCollisionCheck(body1,body2) {

        let dist = Math3D.distance(body1.position,body2.position);

        return dist < (body1.collisionRadius + body2.collisionRadius);
    }

    //Check if point is inside the box volume
    isPointInsideBox(point,body1) {

        //should precompute these for speed with Box objects as reference
        let body1minX = (body1.position[0]-body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1maxX = (body1.position[0]+body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1minY = (body1.position[1]-body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1maxY = (body1.position[1]+body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1minZ = (body1.position[2]-body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1maxZ = (body1.position[2]+body1.collisionRadius)*body1.collisionBoundsScale[0];

        return  (point[0] >= body1minX && point[0] <= body1maxX) &&
                (point[1] >= body1minY && point[1] <= body1maxY) &&
                (point[2] >= body1minZ && point[2] <= body1maxZ);

    }

    //Collision between two axis-aligned boxes. TODO: account for rotation with simple trig modifiers
    boxCollisionCheck(body1,body2) {

        let body1minX = (body1.position[0]-body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1maxX = (body1.position[0]+body1.collisionRadius)*body1.collisionBoundsScale[0];
        let body1minY = (body1.position[1]-body1.collisionRadius)*body1.collisionBoundsScale[1];
        let body1maxY = (body1.position[1]+body1.collisionRadius)*body1.collisionBoundsScale[1];
        let body1minZ = (body1.position[2]-body1.collisionRadius)*body1.collisionBoundsScale[2];
        let body1maxZ = (body1.position[2]+body1.collisionRadius)*body1.collisionBoundsScale[2];

        let body2minX = (body2.position[0]-body2.collisionRadius)*body1.collisionBoundsScale[0];
        let body2maxX = (body2.position[0]+body2.collisionRadius)*body1.collisionBoundsScale[0];
        let body2minY = (body2.position[1]-body2.collisionRadius)*body1.collisionBoundsScale[1];
        let body2maxY = (body2.position[1]+body2.collisionRadius)*body1.collisionBoundsScale[1];
        let body2minZ = (body2.position[2]-body2.collisionRadius)*body1.collisionBoundsScale[2];
        let body2maxZ = (body2.position[2]+body2.collisionRadius)*body1.collisionBoundsScale[2];

        return  (
                    ((body1maxX <= body2maxX && body1maxX >= body2minX) || (body1minX <= body2maxX && body1minX >= body2minX)) &&
                    ((body1maxY <= body2maxY && body1maxY >= body2minY) || (body1minY <= body2maxY && body1minY >= body2minY)) &&
                    ((body1maxZ <= body2maxZ && body1maxZ >= body2minZ) || (body1minZ <= body2maxZ && body1minZ >= body2minZ))
                );
    }

    sphereBoxCollisionCheck(sphere, box) {
        
        let boxMinX = (box.position[0]-box.collisionRadius)*box.collisionBoundsScale[0];
        let boxMaxX = (box.position[0]+box.collisionRadius)*box.collisionBoundsScale[0];
        let boxMinY = (box.position[1]-box.collisionRadius)*box.collisionBoundsScale[1];
        let boxMaxY = (box.position[1]+box.collisionRadius)*box.collisionBoundsScale[1];
        let boxMinZ = (box.position[2]-box.collisionRadius)*box.collisionBoundsScale[2];
        let boxMaxZ = (box.position[2]+box.collisionRadius)*box.collisionBoundsScale[2];

        //let direction = Math.makeVec(sphere.position,box.position);

        //Get closest point to sphere center
        let clamp = [
            Math.max(boxMinX, Math.min(sphere.position[0], boxMaxX)),
            Math.max(boxMinY, Math.min(sphere.position[1], boxMaxY)),
            Math.max(boxMinZ, Math.min(sphere.position[2], boxMaxZ))
        ];

        let dist = Math3D.distance(sphere.position,clamp);

        return dist > sphere.collisionRadius;

    }

    //Closest point on a line from a test point
    closestPointOnLine(testpoint=[0,1,2],point1=[3,4,5],point2=[6,7,8]) {
        let a = [point2[0]-point1[0],point2[1]-point1[1],point2[2]-point1[2]];
        let b = [point1[0]-testpoint[0],point1[1]-testpoint[1],point1[2]-testpoint[2]];
        let c = [point2[0]-testpoint[0],point2[1]-testpoint[1],point2[2]-testpoint[2]];
        let bdota = Math3D.dot(b,a);
        if(bdota <= 0) return point1;
        let cdota = Math3D.dot(c,a);
        if(cdota <= 0) return point2;
        let _bdotapluscdota = 1/(bdota+cdota);
        return [
            point1[0] + ((point2[0]-point1[0])*bdota)*_bdotapluscdota,
            point1[1] + ((point2[1]-point1[1])*bdota)*_bdotapluscdota,
            point1[2] + ((point2[2]-point1[2])*bdota)*_bdotapluscdota
        ];

    }

    //point a,b,c of triangle (following right hand rule for the orientation)
    closestPointOnPolygon(point=[0,1,2], a=3,b=4,c=5) {
        //Find the normal of the polygon
        let n = Math3D.calcNormal(a,b,c);
        //Find the distance from point to the plane given the normal
        let dist = Math3D.dot(point,n) - Math3D.dot(a,n);
        //project p onto the plane by stepping from p to the plane
        let projection = Math3D.vecadd(p,Math3D.vecscale(n,-dist));
        
        //compute edge vectors
        let v0x = c[0] - a[0];
        let v0y = c[1] - a[1];
        let v0z = c[2] - a[2];
        let v1x = b[0] - a[0];
        let v1y = b[1] - a[1];
        let v1z = b[2] - a[2];
        let v2x = projection[0] - a[0];
        let v2y = projection[1] - a[1];
        let v2z = projection[2] - a[2];

        //compute dots
        let dot00 = v0x*v0x+v0y*v0y+v0z*v0z;
        let dot01 = v0x*v1x+v0y*v1y+v0z*v1z;
        let dot02 = v0x*v2x+v0y*v2y+v0z*v2z;
        let dot11 = v1x*v1x+v1y*v1y+v1z*v1z;
        let dot12 = v1x*v2x+v1y*v2y+v1z*v2z;

        //compute barycentric coords (uvs) of projection point
        let denom = dot00*dot11 - dot01*dot01;
        if(Math.abs(denom) < 1e-30) {
            return undefined; //unusable
        }
        let _denom = 1/denom;
        let u = (dot11*dot02 - dot01*dot12)*_denom;
        let v = (dot00*dot12 - dot01*dot02)*_denom;

        //check uv coordinates for validity
        if((u >= 0) && (v >= 0) && (u+v<1)) {
            return projection;
        } else return undefined; //nearest orthogonal point is outside triangle

    }
}