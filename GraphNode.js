import {Math3D} from './Math3D'
//Use this to organize 3D models hierarchically if needed and apply transforms (not very optimal for real time)
export class GraphNode { 
    constructor(parent=null, children=[null], id=null) {
        this.id = id;
        this.parent = parent; //Access/inherit parent object
        this.children = children; //List of child objects for this node, each with independent data
        this.globalPos = {x:0,y:0,z:0}; //Global x,y,z position
        this.localPos = {x:0,y:0,z:0};  //Local x,y,z position offset. Render as global + local pos
        this.globalRot = {x:0,y:0,z:0}; //Global x,y,z rotation (rads)
        this.localRot = {x:0,y:0,z:0}; //Local x,y,z rotation (rads). Render as global + local rot
        this.globalScale = {x:1,y:1,z:1};
        this.localScale = {x:1,y:1,z:1};
        this.functions = []; // List of functions. E.g. function foo(x) {return x;}; this.functions["foo"] = foo; this.functions.foo = foo(x) {return x;}. Got it?

        //3D Rendering stuff
        this.model = null; //
        this.mesh = [0,0,0,1,1,1,1,0,0,0,0,0]; // Model vertex list, array of vec3's xyz, so push x,y,z components. For ThreeJS use THREE.Mesh(vertices, material) to generate a model from this list with the selected material
        this.normals = [];
        this.colors = [0,0,0,255,255,255,255,0,0,0,0,0]; // Vertex color list, array of vec3's rgb or vec4's rgba for outside of ThreeJS. For ThreeJS use THREE.Color("rgb(r,g,b)") for each array item.
        this.materials = []; // Array of materials maps i.e. lighting properties and texture maps.
        this.textures = []; // Array of texture image files.

        if(parent !== null){
            this.inherit(parent);
        }
    }

    inherit(parent) { //Sets globals to be equal to parent info and adds parent functions to this node.
        this.parent = parent;
        this.globalPos = parent.globalPos;
        this.globalRot = parent.globalRot;
        this.globalScale = parent.globalScale;
        this.functions.concat(parent.functions);
        this.children.forEach((child)=>{
            child.inherit(parent);
        });
    }

    addChild(child){ //Add child node reference
        this.children.push(child); //Remember: JS is all pass by object reference.
    }

    removeChild(id){ //Remove child node reference
        this.children.forEach((child, idx)=>{
            if(child.id == id){
                this.children.splice(idx,1);
            }
        });
    }

    translateMeshGlobal(offset=[0,0,0]){ //Recursive global translation of this node and all children
        this.globalPos.x+=offset[0];
        this.globalPos.y+=offset[1];
        this.globalPos.z+=offset[2];

        if(this.mesh.length > 0) this.mesh = Math3D.translateMesh(this.mesh,this.globalPos.x, this.globalPos.y, this.globalPos.z);
        if(this.normals.length > 0) this.normals = Math3D.translateMesh(this.normals,this.globalPos.x, this.globalPos.y, this.globalPos.z);

        this.children.forEach((child)=>{
            child.translateMeshGlobal(offset);
        });
    }

    rotateMeshGlobal(offset=[0,0,0]){ //Offsets the global rotation of this node and all child nodes (radian)
        this.globalRot.x+=offset[0];
        this.globalRot.y+=offset[1];
        this.globalRot.z+=offset[2];

        if(this.mesh.length > 0) this.mesh = Math3D.rotateMesh(this.mesh,this.globalRot.x, this.globalRot.y, this.globalRot.z);
        if(this.normals.length > 0) this.normals = Math3D.rotateMesh(this.normals,this.globalRot.x, this.globalRot.y, this.globalRot.z);

        this.children.forEach((child)=>{
            child.rotateMeshGlobal(offset);
        });
    }

    translateMeshLocal(offset=[0,0,0]){ //Recursive global translation of this node and all children
        this.localPos.x+=offset[0];
        this.localPos.y+=offset[1];
        this.localPos.z+=offset[2];

        if(this.mesh.length > 0) this.mesh = Math3D.translateMesh(this.mesh,this.localPos.x, this.localPos.y, this.localPos.z);
        if(this.normals.length > 0) this.normals = Math3D.translateMesh(this.normals,this.localPos.x, this.localPos.y, this.localPos.z);

        this.children.forEach((child)=>{
            child.translateMeshLocal(offset);
        });
    }

    translateMeshGlobal(offset=[0,0,0]){ //Offsets the global rotation of this node and all child nodes (radian)
        this.localRot.x+=offset[0];
        this.localRot.y+=offset[1];
        this.localRot.z+=offset[2];

        if(this.mesh.length > 0) this.mesh = Math3D.rotateMesh(this.mesh,this.globalPos.x, this.globalPos.y, this.globalPos.z);
        if(this.normals.length > 0) this.normals = Math3D.rotateMesh(this.normals,this.globalPos.x, this.globalPos.y, this.globalPos.z);

        this.children.forEach((child)=>{
            child.translateMeshGlobal(offset);
        });
    }

    scaleMeshLocal(scalar=[1,1,1]){
        this.localScale.x+=scalar[0];
        this.localScale.y+=scalar[1];
        this.localScale.z+=scalar[2];

        if(this.mesh.length > 0) this.mesh = Math3D.scaleMesh(this.mesh,this.localScale.x, this.localScale.y, this.localScale.z);
        if(this.normals.length > 0) this.normals = Math3D.scaleMesh(this.normals,this.localScale.x, this.localScale.y, this.localScale.z);

        this.children.forEach((child)=>{
            child.scaleMeshLocal(offset);
        });
    }

    scaleMeshGlobal(scalar=[1,1,1]){
        this.globalScale.x+=scalar[0];
        this.globalScale.y+=scalar[1];
        this.globalScale.z+=scalar[2];

        if(this.mesh.length > 0) this.mesh = Math3D.scaleMesh(this.mesh,this.globalScale.x, this.globalScale.y, this.globalScale.z);
        if(this.normals.length > 0) this.normals = Math3D.scaleMesh(this.normals,this.globalScale.x, this.globalScale.y, this.globalScale.z);

        this.children.forEach((child)=>{
            child.scaleMeshGlobal(offset);
        });
    }

    applyMeshTransforms() { //Get mesh with rotation and translation applied
        var rotated = Math3D.rotateMesh(this.mesh,this.globalRot.x+this.localRot.x,this.globalRot.y+this.localRot.y,this.globalRot.z+this.localRot.z);
        var translated = Math3D.translateMesh(rotated,this.globalPos.x+this.localPos.x, this.globalPos.y+this.localPos.y, this.globalPos.z+this.localPos.z);
        var scaled = Math3D.scaleMesh(translated, this.globalScale.x+this.localScale.x, this.globalScale.y+this.localScale.y, this.globalScale.z+this.localScale.z);
        return scaled;
    }

}
