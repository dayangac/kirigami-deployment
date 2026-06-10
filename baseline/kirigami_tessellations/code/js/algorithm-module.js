// js/algorithm-module.js
export class AlgorithmModule {
    constructor() {
        this.meshData = null;
        this.patternData = null;
        this.initialized = false;
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.mesh = null;
        this.wireframe = null;
        this.showWireframe = false;
        this.currentGroundFaces = []
        this.showDeployedWireframe = false;
    }
    

    
    toggleWireframe() {
        this.showWireframe = !this.showWireframe;
        
        const toggleWireframeBtn = document.getElementById('toggleWireframeBtn');
        if (toggleWireframeBtn) {
            // Use fixed text length with spaces to maintain same width
            toggleWireframeBtn.textContent = this.showWireframe ? 
                'show mesh' : 'show wireframe';
            
            // Add/remove active class for visual feedback
            if (this.showWireframe) {
                toggleWireframeBtn.classList.add('active');
            } else {
                toggleWireframeBtn.classList.remove('active');
            }
        }
        
        this.updateWireframeVisibility();
    }

    
    updateWireframeVisibility() {
        if (!this.mesh || !this.wireframe) return;
        
        this.wireframe.visible = this.showWireframe;
        
        // Also update the main mesh material to be more transparent when wireframe is shown
        if (this.mesh.material) {
            this.mesh.material.wireframe = this.showWireframe;
            this.mesh.material.opacity = this.showWireframe ? 0.3 : 1.0;
            this.mesh.material.transparent = this.showWireframe;
            this.mesh.material.needsUpdate = true;
        }
    }
    
    initDeployedConfigThreeJS() {
      const container = document.getElementById('canvasOptiDeploy');
      // Clear previous content
      container.innerHTML = '';
      this.deployedScene = new THREE.Scene();
      this.deployedScene.background = new THREE.Color(0xf9f9f9);
      this.deployedScene.add(new THREE.AmbientLight(0x606060));
      const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
      directionalLight.position.set(1, 1, 1).normalize();
      this.deployedScene.add(directionalLight);
      this.deployedCamera = new THREE.PerspectiveCamera(75, container.clientWidth / container.clientHeight, 0.1, 1000);
      this.deployedCamera.position.z = 5;
      this.deployedRenderer = new THREE.WebGLRenderer({ antialias: true });
      this.deployedRenderer.setSize(container.clientWidth, container.clientHeight);
      this.deployedRenderer.setPixelRatio(window.devicePixelRatio);
      container.appendChild(this.deployedRenderer.domElement);
      this.deployedControls = new THREE.OrbitControls(this.deployedCamera, this.deployedRenderer.domElement);
      this.deployedControls.enableDamping = true;
      this.deployedControls.dampingFactor = 0.25;
    }
    initThreeJS() {
        const container = document.getElementById('meshViewer');
        if (!container) return;
        
        // Clear previous content
        container.innerHTML = '';
        
        // Create Three.js scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0xf0f0f0);
        
        // Add some basic lighting
        const ambientLight = new THREE.AmbientLight(0x606060);
        this.scene.add(ambientLight);
        
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(1, 1, 1).normalize();
        this.scene.add(directionalLight);
        
        // Create camera
        this.camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);
        this.camera.position.z = 5;
        
        // Create renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        container.appendChild(this.renderer.domElement);
        
        // Add orbit controls for interactive camera
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.25;
        
        // Add grid helper
        const gridHelper = new THREE.GridHelper(10, 10);
        this.scene.add(gridHelper);
        
        // Add axes helper
        const axesHelper = new THREE.AxesHelper(5);
        this.scene.add(axesHelper);
        
        // Handle window resize
        window.addEventListener('resize', () => this.onWindowResize());
        
        // Start animation loop
        this.animate();
    }
    
    onWindowResize() {
        const container = document.getElementById('meshViewer');
        if (!container || !this.camera || !this.renderer) return;
        
        this.camera.aspect = container.clientWidth / container.clientHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(container.clientWidth, container.clientHeight);
    }
    
    animate() {
        requestAnimationFrame(() => this.animate());
        
        if (this.controls) {
            this.controls.update();
        }
        
        if (this.renderer && this.scene && this.camera) {
            this.renderer.render(this.scene, this.camera);
        }
        if (this.deployedRenderer && this.deployedScene && this.deployedCamera) {
          this.deployedRenderer.render(this.deployedScene, this.deployedCamera);
        }
    }
    
    handleMeshFileSelect(event) {
        const file = event.target.files[0];
        if (!file) return;
        
        const reader = new FileReader();
        const extension = file.name.split('.').pop().toLowerCase();
        
        reader.onload = (e) => {
            try {
                if (extension === 'obj') {
                    this.meshData = this.parseOBJ(e.target.result);
                    document.getElementById('initOptimizationButton').disabled = false;
                } else if (extension === 'off') {
                    this.meshData = this.parseOFF(e.target.result);
                    document.getElementById('initOptimizationButton').disabled = false;
                } else {
                    throw new Error('Unsupported file format');
                }
                
                this.displayMesh(this.meshData);
                this.updatePrintInfo();
                
            } catch (error) {
                console.error('Error loading mesh:', error);
                alert('Error loading mesh: ' + error.message);
            }
        };
        
        reader.readAsText(file);
    }
    
    parseOBJ(objString) {
        const vertices = [];
        const faces = [];
        
        const lines = objString.split('\n');
        lines.forEach(line => {
            const parts = line.trim().split(/\s+/);
            if (parts.length === 0) return;
            
            const type = parts[0];
            
            if (type === 'v') {
                // Vertex
                vertices.push([
                    parseFloat(parts[1]),
                    parseFloat(parts[2]),
                    parseFloat(parts[3])
                ]);
            } else if (type === 'f') {
                // Face
                const faceVertices = [];
                for (let i = 1; i < parts.length; i++) {
                    const vertexIndex = parseInt(parts[i].split('/')[0]) - 1;
                    if (!isNaN(vertexIndex)) {
                        faceVertices.push(vertexIndex);
                    }
                }
                if (faceVertices.length >= 3) {
                    faces.push(faceVertices);
                }
            }
        });
        
        return { vertices, faces };
    }
    
    parseOFF(offString) {
        const lines = offString.split('\n').filter(line => 
            line.trim() !== '' && !line.startsWith('#')
        );
        
        if (lines[0].trim() !== 'OFF') {
            throw new Error('Invalid OFF file format');
        }
        
        const header = lines[1].trim().split(/\s+/);
        const vertexCount = parseInt(header[0]);
        const faceCount = parseInt(header[1]);
        
        const vertices = [];
        const faces = [];
        
        // Parse vertices
        for (let i = 2; i < 2 + vertexCount; i++) {
            const coords = lines[i].trim().split(/\s+/);
            vertices.push([
                parseFloat(coords[0]),
                parseFloat(coords[1]),
                parseFloat(coords[2])
            ]);
        }
        
        // Parse faces
        for (let i = 2 + vertexCount; i < 2 + vertexCount + faceCount; i++) {
            const faceData = lines[i].trim().split(/\s+/);
            const vertexCountInFace = parseInt(faceData[0]);
            const faceVertices = [];
            
            for (let j = 1; j <= vertexCountInFace; j++) {
                faceVertices.push(parseInt(faceData[j]));
            }
            
            if (vertexCountInFace >= 3) {
                faces.push(faceVertices);
            }
        }
        
        return { vertices, faces };
    }
    
    displayMesh(meshData) {
        if (!this.scene || !meshData) return;
        
        // Remove existing mesh and wireframe if any
        if (this.mesh) {
            this.scene.remove(this.mesh);
        }
        if (this.wireframe) {
            this.scene.remove(this.wireframe);
        }
        
        // Create geometry from parsed data
        const geometry = new THREE.BufferGeometry();
        
        // Create vertices array
        const vertices = [];
        meshData.vertices.forEach(vertex => {
            vertices.push(vertex[0], vertex[1], vertex[2]);
        });
        
        // Create faces/indices array
        const indices = [];
        meshData.faces.forEach(face => {
            // Convert to triangles (assuming convex polygons)
            for (let i = 1; i < face.length - 1; i++) {
                indices.push(face[0], face[i], face[i + 1]);
            }
        });
        
        // Set geometry attributes
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setIndex(indices);
        geometry.computeVertexNormals();
        
        // Create material and mesh
        const material = new THREE.MeshPhongMaterial({
            color: 0x78c6a3, // Theme color
            side: THREE.DoubleSide,
            flatShading: true,
            wireframe: false
        });
        
        this.mesh = new THREE.Mesh(geometry, material);
        this.scene.add(this.mesh);
        this.updatePrintInfo();

        // Create wireframe
        const wireframeMaterial = new THREE.LineBasicMaterial({
            color: 0x000000,
            linewidth: 5,
            depthTest: false // Ensure wireframe is always on top
        });
        
        const wireframeGeometry = new THREE.WireframeGeometry(geometry);
        this.wireframe = new THREE.LineSegments(wireframeGeometry, wireframeMaterial);
        this.wireframe.visible = this.showWireframe;
        this.scene.add(this.wireframe);
        
        // Center and scale the mesh
        geometry.computeBoundingBox();
        const center = new THREE.Vector3();
        geometry.boundingBox.getCenter(center);
        this.mesh.position.sub(center);
        this.wireframe.position.sub(center);
        
        // Auto-rotate for better viewing
        if (this.controls) {
            this.controls.autoRotate = true;
            this.controls.autoRotateSpeed = 1.0;
        }
        
        // Fit camera to mesh
        this.fitCameraToObject(this.mesh, this.camera, this.controls);
        
        // Update wireframe button state
        const toggleWireframeBtn = document.getElementById('toggleWireframeBtn');
        if (toggleWireframeBtn) {
            toggleWireframeBtn.disabled = false;
            // Use fixed text length with spaces to maintain same width
            toggleWireframeBtn.textContent = this.showWireframe ? 
                'show mesh' : 'show wireframe';
            
            // Add/remove active class for visual feedback
            if (this.showWireframe) {
                toggleWireframeBtn.classList.add('active');
            } else {
                toggleWireframeBtn.classList.remove('active');
            }
        }
    }
    
    fitCameraToObject(object, camera, controls) {
      const boundingBox = new THREE.Box3().setFromObject(object);
      const center = new THREE.Vector3();
      boundingBox.getCenter(center);

      const size = boundingBox.getSize(new THREE.Vector3());
      const maxDim = Math.max(size.x, size.y, size.z);
      const fov = camera.fov * (Math.PI / 180);

      let distance = Math.abs(maxDim / (2 * Math.tan(fov / 2)));
      distance *= 1.5; // padding

      // figure out camera’s current direction
      const dir = new THREE.Vector3();
      camera.getWorldDirection(dir);
      dir.normalize();

      // place camera along that direction, away from center
      camera.position.copy(center.clone().add(dir.clone().multiplyScalar(-distance)));

      // update controls target
      controls.target.copy(center);
      controls.update();
  }

    addFaceCountToCanvas(canvas, faceCount) {
      const ctx = canvas.getContext('2d');
      
      // Clear previous text by redrawing the pattern
      window.previewModule.drawPatternOnCanvas(
          { verts: this.currentGroundVerts || [], faces: this.currentGroundFaces || [] }, 
          false, 
          canvas
      );
      
      // Draw face count at the bottom center, matching the mesh stats style
      ctx.fillStyle = '#444';
      ctx.font = '40px Arial';
      ctx.textAlign = 'right';
      ctx.textBaseline = 'top';
      ctx.fillText(`#faces: ${faceCount}`, canvas.width -10 , 10);
  }
        
    updatePrintInfo() {
        // Update the stats panel below the viewer
        const statsPanel = document.getElementById('meshStats'); // stats container below mesh viewers

        if (!statsPanel || !this.meshData) return;

        if (statsPanel) {
            statsPanel.innerHTML = `
                #vert: <strong>${this.meshData.vertices.length}</strong>, #face:<strong>${this.meshData.faces.length}</strong>
            `;
        }
    }

    clearStats() {
        const statsPanel = document.getElementById('meshStats');
        if (statsPanel) {
            statsPanel.innerHTML = 'no mesh loaded';
        }
        
        const infoPanel = document.getElementById('printInfo');
        if (infoPanel) {
            infoPanel.innerHTML = '<p>Load a mesh to see print information.</p>';
        }
    }
    
    initThreeJS(container) {
        if (!container) return;
        
        // Clear previous content
        container.innerHTML = '';
        
        // Create Three.js scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0xf0f0f0);
        
        // Add some basic lighting
        const ambientLight = new THREE.AmbientLight(0x606060);
        this.scene.add(ambientLight);
        
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(1, 1, 1).normalize();
        this.scene.add(directionalLight);
        
        // Create camera
        this.camera = new THREE.PerspectiveCamera(75, container.clientWidth / container.clientHeight, 0.1, 1000);
        this.camera.position.z = 5;
        
        // Create renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        container.appendChild(this.renderer.domElement);
        
        // Add orbit controls for interactive camera
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.25;
        
        // Add grid helper
        const gridHelper = new THREE.GridHelper(10, 10);
        this.scene.add(gridHelper);
        
        // Add axes helper
        const axesHelper = new THREE.AxesHelper(5);
        this.scene.add(axesHelper);
        
        // Handle window resize
        window.addEventListener('resize', () => this.onWindowResize());
        
        // Start animation loop
        this.animate();
    }

    updateGroundConfigVisualization(verts, faces) {
        // Store current ground vertices and faces for display
        this.currentGroundVerts = verts;
        this.currentGroundFaces = faces;
      
        let pattern = { verts: verts, faces: faces };
        window.previewModule.drawPatternOnCanvas(pattern, false, this.initConfigCanvas);
        
        // Add face count display to the canvas
        this.addFaceCountToCanvas(this.initConfigCanvas, faces.length);
    }

    // updateLiftedConfigVisualization(verts, faces) {
    //   if (this.deployedMesh) {
    //     this.deployedScene.remove(this.deployedMesh);
    //     this.deployedMesh.geometry.dispose();
    //     this.deployedMesh.material.dispose();
    //     this.deployedMesh = null;
    //   }
    //   if (verts.length === 0 || faces.length === 0)
    //     return;
    //   const geometry = new THREE.BufferGeometry();
    //   const vertices = [];
    //   verts.forEach(v => {
    //     vertices.push(v[0], v[1], v[2]);
    //   });
    //   const indices = [];
    //   faces.forEach(face => {
    //     for (let i = 1; i < face.length - 1; i++) {
    //       indices.push(face[0], face[i], face[i + 1]);
    //     }
    //   });
    //   geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
    //   geometry.setIndex(indices);
    //   geometry.computeVertexNormals();
    //   const material = new THREE.MeshPhongMaterial({
    //     color: 0x78c6a3,
    //     side: THREE.DoubleSide,
    //     flatShading: true,
    //     wireframe: false
    //   });
    //   this.deployedMesh = new THREE.Mesh(geometry, material);
    //   this.deployedScene.add(this.deployedMesh);
    //   geometry.computeBoundingBox();
    //   const center = new THREE.Vector3();
    //   geometry.boundingBox.getCenter(center);
    //   this.deployedMesh.position.sub(center);
    //   if (this.deployedControls) {
    //     this.deployedControls.autoRotate = true;
    //     this.deployedControls.autoRotateSpeed = 1.0;
    //   }
    //   this.fitCameraToObject(this.deployedMesh, this.deployedCamera, this.deployedControls);
    // }

    updateLiftedConfigVisualization(verts, faces) {
      if (this.deployedMesh) {
          this.deployedScene.remove(this.deployedMesh);
          this.deployedMesh.geometry.dispose();
          this.deployedMesh.material.dispose();
          this.deployedMesh = null;
      }
      if (this.deployedWireframe) {
          this.deployedScene.remove(this.deployedWireframe);
          this.deployedWireframe.geometry.dispose();
          this.deployedWireframe.material.dispose();
          this.deployedWireframe = null;
      }
      
      if (verts.length === 0 || faces.length === 0)
          return;
          
      // Create the main mesh (triangulated for rendering)
      const geometry = new THREE.BufferGeometry();
      const vertices = [];
      verts.forEach(v => {
          vertices.push(v[0], v[1], v[2]);
      });
      
      const indices = [];
      faces.forEach(face => {
          for (let i = 1; i < face.length - 1; i++) {
              indices.push(face[0], face[i], face[i + 1]);
          }
      });
      
      geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
      geometry.setIndex(indices);
      geometry.computeVertexNormals();
      
      const material = new THREE.MeshPhongMaterial({
          color: 0x78c6a3,
          side: THREE.DoubleSide,
          flatShading: true,
          wireframe: false
      });
      
      this.deployedMesh = new THREE.Mesh(geometry, material);
      this.deployedScene.add(this.deployedMesh);
      
      // Create wireframe that's always visible on top
      const edgeGeometry = new THREE.BufferGeometry();
      const edgeVertices = [];
      
      faces.forEach(face => {
          const vertexCount = face.length;
          for (let i = 0; i < vertexCount; i++) {
              const v1 = verts[face[i]];
              const v2 = verts[face[(i + 1) % vertexCount]];
              
              edgeVertices.push(v1[0], v1[1], v1[2]);
              edgeVertices.push(v2[0], v2[1], v2[2]);
          }
      });
      
      edgeGeometry.setAttribute('position', new THREE.Float32BufferAttribute(edgeVertices, 3));
      
      const wireframeMaterial = new THREE.LineBasicMaterial({
          color: 0x000000,
          linewidth: 20,
      });
      
      this.deployedWireframe = new THREE.LineSegments(edgeGeometry, wireframeMaterial);
      this.deployedScene.add(this.deployedWireframe); // Always add to scene
      
      geometry.computeBoundingBox();
      const center = new THREE.Vector3();
      geometry.boundingBox.getCenter(center);
      this.deployedMesh.position.sub(center);
      this.deployedWireframe.position.sub(center);
      
      if (this.deployedControls) {
          this.deployedControls.autoRotate = true;
          this.deployedControls.autoRotateSpeed = 1.0;
      }
      this.fitCameraToObject(this.deployedMesh, this.deployedCamera, this.deployedControls);
  }
    
    

    init() {
        if (this.initialized) return;
        this.initialized = true;
        
        const meshViewer = document.getElementById('meshViewer');
        // Initialize Three.js with provided container
        this.initDeployedConfigThreeJS();
        this.initThreeJS(meshViewer);

        // Clear stats on initialization
        this.clearStats();
        
        console.log('Algorithm module initialized with Three.js');

        // Set up pattern choosing UI.
        this.unitPatternCanvas = document.getElementById('previewCanvasP4');
        const dpr = window.devicePixelRatio || 1;
        this.unitPatternCanvas.width = this.unitPatternCanvas.clientWidth * dpr * 2;
        this.unitPatternCanvas.height = this.unitPatternCanvas.clientHeight * dpr * 2;
        this.initConfigCanvas = document.getElementById('initialConfigCanvas');
        this.initConfigCanvas.width = this.initConfigCanvas.clientWidth * dpr * 2;
        this.initConfigCanvas.height = this.initConfigCanvas.clientHeight * dpr * 2;
        document.getElementById('patternPresetSelectP4').addEventListener('change', async (e) => {
          let groups = await window.designModule.loadPresetPatternGroups(e.target.value);
          await window.kiriCallAsync('set_pattern', groups);
          let cut_pattern = await window.kiriCallAsync('open_pattern', [Math.PI / 10, 3, 3]);
          window.previewModule.drawPatternOnCanvas(cut_pattern, false, this.unitPatternCanvas);
        });

        // Update errors.
        let update_errors = (result) => {
          // Check if maxError is checked.
          let maxError = document.getElementById('maxError').checked;
          let planar_err, shape_err, rigid_err;
          if (maxError) {
            planar_err = result['planarity_max'];
            shape_err = result['close_max']
            rigid_err = result['rigid_max']
          } else {
            planar_err = result['planarity_avg'];
            shape_err = result['close_avg']
            rigid_err = result['rigid_avg']
          }
          let planar_err_span = document.getElementById('planarityError');
          let shape_err_span = document.getElementById('shapeError');
          let rigid_err_span = document.getElementById('rigidityError');
          // Update text content with 4 decimal places.
          planar_err_span.textContent = planar_err.toFixed(4) + "°";
          if (planar_err > 0.5) {
            planar_err_span.style.color = "red";
          } else {
            planar_err_span.style.color = "";
          }
          shape_err_span.textContent = shape_err.toFixed(4);
          rigid_err_span.textContent = (rigid_err * 100).toFixed(4) + "%";
          if (rigid_err > 0.01) {
            rigid_err_span.style.color = "red";
          } else {
            rigid_err_span.style.color = "";
          }
        };
        document.getElementById('maxError').addEventListener('change', (e) => {
          window.kiriCallAsync('get_errors', []).then(update_errors);
        });
        document.getElementById('avgError').addEventListener('change', (e) => {
          window.kiriCallAsync('get_errors', []).then(update_errors);
        });

        // Optimization UI.
        let optimization_initialized = false;
        document.getElementById('initOptimizationButton').addEventListener('click', async () => {
          window.kiriCallAsync('init_opt', [window.algorithmModule.meshData.vertices, window.algorithmModule.meshData.faces]).then(result => {
            this.updateGroundConfigVisualization(result['ground_verts'], result['ground_faces']);
            this.updateLiftedConfigVisualization(result['lifted_verts'], result['lifted_faces']);
            document.getElementById('runOptimizationBtn').disabled = false;
            optimization_initialized = true;
            document.getElementById('scaleInput').value = 1;
            document.getElementById('rotationInput').value = 0;
            window.kiriCallAsync('get_errors', []).then(update_errors);
          });
        });

        // Lift params UI.
        let already_calculating = false;
        let update_params_func = () => {
          if (already_calculating) return;
          already_calculating = true;
          // Get params.
          let scale = parseFloat(document.getElementById('scaleInput').value);
          let angle = parseFloat(document.getElementById('rotationInput').value) * Math.PI / 180;
          window.kiriCallAsync('update_lift_params', [scale, angle]).then(result => {
            this.updateGroundConfigVisualization(result['ground_verts'], result['ground_faces']);
            this.updateLiftedConfigVisualization(result['lifted_verts'], result['lifted_faces']);

            // Check if need to calculate again.
            let new_scale = parseFloat(document.getElementById('scaleInput').value);
            let new_angle = parseFloat(document.getElementById('rotationInput').value) * Math.PI / 180;
            already_calculating = false;
            if (new_scale !== scale || new_angle !== angle) {
              setTimeout(update_params_func, 10);
            } else {
              window.kiriCallAsync('get_errors', []).then(update_errors);
            }
          });
        };
        document.getElementById('scaleInput').addEventListener('input', (e) => {
          if (!optimization_initialized || already_calculating) return;
          update_params_func();
        });
        document.getElementById('rotationInput').addEventListener('input', (e) => {
          if (!optimization_initialized || already_calculating) return;
          update_params_func();
        });

        // Run optimization button.
        let running_opt = false;
        let run_opt = () => {
          window.kiriCallAsync('optimize', []).then(result => {
            this.updateGroundConfigVisualization(result['ground_verts'], result['ground_faces']);
            this.updateLiftedConfigVisualization(result['lifted_verts'], result['lifted_faces']);
            // Update errors.
            update_errors(result);
            // Continue optimization if still running.
            if (running_opt) {
              setTimeout(run_opt, 10);
            }
          });
        };
        document.getElementById('runOptimizationBtn').addEventListener('click', async () => {
          running_opt = !running_opt;
          if (running_opt) {
            run_opt();
            document.getElementById('runOptimizationBtn').textContent = 'Stop Optimization';
            document.getElementById('runOptimizationBtn').style.backgroundColor = "#f44336"; // Red color
          } else {
            document.getElementById('runOptimizationBtn').textContent = 'Run Optimization';
            document.getElementById('runOptimizationBtn').style.backgroundColor = ""; // Reset to default
          }
        });

    }


}

// Initialize when DOM is loaded
// document.addEventListener('DOMContentLoaded', () => {
//     window.algorithmModule = new AlgorithmModule();
//     window.algorithmModule.init();
// });