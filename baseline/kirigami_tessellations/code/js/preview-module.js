// js/preview-module.js
export class PreviewModule {
    constructor(canvasId, statusId, designModule) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');
        this.statusElement = document.getElementById(statusId);
        this.designModule = designModule;
        
        this.basePattern = null;
        this.cutOpenPattern = null;
        this.unitPattern = null;
        this.isCutOpen = false;
        
        // Zoom and pan state
        this.zoom = 1;
        this.panX = 0;
        this.panY = 0;
        this.isPanning = false;
        this.lastPanX = 0;
        this.lastPanY = 0;


        
        // Animation state
        this.animationId = null;
        this.isAnimating = false;
        this.animationDirection = 1; // 1 for increasing, -1 for decreasing
        this.animationSpeed = 1; // degrees per frame

        // Add deployment state
        this.maxAreaAngle = 0;
        this.maxOpeningAngle = 0;
        
        
        this.setupEventListeners();
        this.updateCanvasDimensions();
    }

    setupEventListeners() {
        // Wheel event for zooming
        this.canvas.addEventListener('wheel', (e) => {
            if (e.shiftKey) {
                e.preventDefault();
                this.handleZoom(e.deltaY, e.offsetX, e.offsetY);
            }
        });

        // Mouse events for panning
        this.canvas.addEventListener('mousedown', (e) => {
            if (e.shiftKey) {
                this.isPanning = true;
                this.lastPanX = e.clientX;
                this.lastPanY = e.clientY;
                this.canvas.style.cursor = 'grabbing';
            }
        });

        this.canvas.addEventListener('mousemove', (e) => {
            if (this.isPanning) {
                const deltaX = e.clientX - this.lastPanX;
                const deltaY = e.clientY - this.lastPanY;
                this.panX += deltaX;
                this.panY += deltaY;
                this.lastPanX = e.clientX;
                this.lastPanY = e.clientY;
                this.drawCurrentPattern();
            }
        });

        this.canvas.addEventListener('mouseup', () => {
            if (this.isPanning) {
                this.isPanning = false;
                this.canvas.style.cursor = 'default';
            }
        });

        this.canvas.addEventListener('mouseleave', () => {
            if (this.isPanning) {
                this.isPanning = false;
                this.canvas.style.cursor = 'default';
            }
        });

        // Opening angle slider
        const openingAngleSlider = document.getElementById('openingAngleSlider');
        const angleValue = document.getElementById('angleValue');
        
        if (openingAngleSlider && angleValue) {
            openingAngleSlider.addEventListener('input', (e) => {
                const value = parseInt(e.target.value);
                angleValue.textContent = `${value}°`;
                
                if (this.cutOpenPattern) {
                    const openingAngleRad = value * Math.PI / 180;
                    this.updatePatternWithOpeningAngle(openingAngleRad);
                }
            });
        }

        // Auto-play toggle
        if (autoPlayToggle) {
            autoPlayToggle.addEventListener('change', (e) => {
                if (e.target.checked) {
                    this.startAnimation();
                } else {
                    this.stopAnimation();
                }
            });
        }

        // Deployment radio buttons
        const deployIni = document.getElementById('deploy-ini');
        const deployMaxArea = document.getElementById('deploy-max-area');
        const deployMaxAngle = document.getElementById('deploy-max-angle');
        
        if (deployIni && deployMaxArea && deployMaxAngle) {
            deployIni.addEventListener('change', (e) => {
                if (e.target.checked) {
                    this.setDeploymentAngle('ini');
                }
            });
            
            deployMaxArea.addEventListener('change', (e) => {
                if (e.target.checked) {
                    this.setDeploymentAngle('max-area');
                }
            });
            
            deployMaxAngle.addEventListener('change', (e) => {
                if (e.target.checked) {
                    this.setDeploymentAngle('max-angle');
                }
            });
        }
        
        
        

        const generatePatternBtn = document.getElementById('generatePatternBtn');
        const angleControl = document.querySelector('.angle-control');
        // const angleControl = document.querySelector('.angle-control');
        
        if (generatePatternBtn) {
            generatePatternBtn.addEventListener('click', () => {
                // Generate the base pattern first
                this.generateBasePattern();
                
                // Show the angle control
                angleControl.style.display = 'flex';
                
                // Show the deployed pattern
                this.cut_open_pattern();

                // Start auto-play animation if enabled
                if (autoPlayToggle.checked) {
                    this.startAnimation();
                }
            });
        }

        // Toggle deployment button
        // const toggleDeploymentBtn = document.getElementById('toggleDeploymentBtn');
        // const angleControl = document.querySelector('.angle-control');
        
        // if (toggleDeploymentBtn && angleControl) {
        //     toggleDeploymentBtn.addEventListener('click', () => {
        //         const isDeployed = angleControl.style.display !== 'none';
                
        //         if (!isDeployed && this.basePattern) {
        //             // Show deployed pattern
        //             this.cut_open_pattern();
        //             angleControl.style.display = 'block';
        //             toggleDeploymentBtn.textContent = 'Show Flat Pattern';
        //         } else {
        //             // Show flat pattern
        //             this.show_base_pattern();
        //             angleControl.style.display = 'none';
        //             toggleDeploymentBtn.textContent = 'Show Deployed Pattern';
        //         }
        //     });
        // }
    }

    setDeploymentAngle(type) {
        const openingAngleSlider = document.getElementById('openingAngleSlider');
        const angleValue = document.getElementById('angleValue');
        const autoPlayToggle = document.getElementById('autoPlayToggle');
        
        if (!openingAngleSlider || !angleValue || !autoPlayToggle) return;
        
        // Turn off auto-play when selecting a deployment option
        autoPlayToggle.checked = false;
        this.stopAnimation();
        
        let angle = 0;
        
        switch (type) {
            case 'ini':
                angle = 0;
                break;
            case 'max-area':
                angle = this.maxAreaAngle || parseInt(openingAngleSlider.max);
                break;
            case 'max-angle':
                angle = this.maxOpeningAngle || parseInt(openingAngleSlider.max);
                break;
        }
        
        // Update slider and pattern
        openingAngleSlider.value = angle;
        angleValue.textContent = `${angle}°`;
        
        const openingAngleRad = angle * Math.PI / 180;
        this.updatePatternWithOpeningAngle(openingAngleRad);
    }

    // Animation methods
    startAnimation() {
        if (this.isAnimating) return;
        
        this.isAnimating = true;
        this.animateOpeningAngle();
    }
    
    stopAnimation() {
        this.isAnimating = false;
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
            this.animationId = null;
        }
    }
    
    animateOpeningAngle() {
        if (!this.isAnimating) return;
        
        const openingAngleSlider = document.getElementById('openingAngleSlider');
        const angleValue = document.getElementById('angleValue');
        
        if (!openingAngleSlider || !this.cutOpenPattern) {
            this.stopAnimation();
            return;
        }
        
        let currentValue = parseInt(openingAngleSlider.value);
        const maxAngle = parseInt(openingAngleSlider.max);
        
        // Update value based on direction
        currentValue += this.animationDirection * this.animationSpeed;
        
        // Change direction if we hit boundaries
        if (currentValue >= maxAngle) {
            currentValue = maxAngle;
            this.animationDirection = -1;
        } else if (currentValue <= 0) {
            currentValue = 0;
            this.animationDirection = 1;
        }
        
        // Update slider and pattern
        openingAngleSlider.value = currentValue;
        angleValue.textContent = `${currentValue}°`;
        
        const openingAngleRad = currentValue * Math.PI / 180;
        this.updatePatternWithOpeningAngle(openingAngleRad);
        
        // Continue animation
        this.animationId = requestAnimationFrame(() => this.animateOpeningAngle());
    }


    updateCanvasDimensions() {
        const container = this.canvas.parentElement;
        const size = Math.min(container.clientWidth, 400);
        this.canvas.width = size;
        this.canvas.height = size;
    }

    generateBasePattern() {
        if (!window.kiri || !this.designModule.groups || this.designModule.groups.length < 2) {
            this.statusElement.textContent = 'Need at least 2 groups to generate pattern';
            return;
        }

        try {
            this.basePattern = window.kiri.createBasePattern(this.designModule.groups);
            if (this.basePattern.error) {
                throw new Error(this.basePattern.error);
            }
            
            this.show_base_pattern();
            this.statusElement.textContent = 'Base pattern generated. Click "Show Deployed Pattern" to see deployment.';
            
        } catch (error) {
            this.statusElement.textContent = `Error: ${error.message}`;
            console.error('Pattern generation error:', error);
        }
    }

    show_base_pattern() {
        if (!this.basePattern) return;
        
        this.isCutOpen = false;
        const repeat_x = 3; // Default repeat values
        const repeat_y = 3;
        
        this.unitPattern = window.kiri.createFinalPattern(this.basePattern, 0, 1, 1);
        this.cutOpenPattern = window.kiri.createFinalPattern(this.basePattern, 0, repeat_x, repeat_y, false);
        let opening_angle = this.unitPattern.max_opening_angle || 0;
        this.maxOpeningAngle = Math.floor(opening_angle * 180 / Math.PI);
        
        // todo: compute the max area angle
        this.maxAreaAngle = this.maxOpeningAngle;
        document.getElementById('openingAngleSlider').max = this.maxOpeningAngle;
        document.getElementById('openingAngleSlider').value = 0;
      

        this.drawPattern(this.cutOpenPattern, false);
        this.statusElement.textContent = 'Flat pattern displayed';

        

        // Hide the angle control when showing base pattern
        const angleControl = document.querySelector('.angle-control');
        if (angleControl) {
            angleControl.style.display = 'flex';
        }

        // Reset deployment options to "ini"
        const deployIni = document.getElementById('deploy-ini');
        if (deployIni) {
            deployIni.checked = true;
        }
    }

    cut_open_pattern() {
        if (!this.basePattern) return;
        
        this.isCutOpen = true;
        const openingAngle = parseInt(document.getElementById('openingAngleSlider').value || 0) * Math.PI / 180;
        const repeat_x = parseInt(document.getElementById('repeatX').value) || 3;
        const repeat_y = parseInt(document.getElementById('repeatY').value) || 3;
        
        this.unitPattern = window.kiri.createFinalPattern(this.basePattern, openingAngle, 1, 1);
        this.cutOpenPattern = window.kiri.createFinalPattern(this.basePattern, openingAngle, repeat_x, repeat_y);  

        this.drawPattern(this.cutOpenPattern, true);
        this.statusElement.textContent = `Deployed pattern at ${(openingAngle * 180/Math.PI).toFixed(1)}°`;
    }

    updatePatternWithOpeningAngle(openingAngleRad) {
        if (!this.basePattern) return;
        
        const repeat_x = parseInt(document.getElementById('repeatX').value) || 3;
        const repeat_y = parseInt(document.getElementById('repeatY').value) || 3;
        
        this.cutOpenPattern = window.kiri.createFinalPattern(this.basePattern, openingAngleRad, repeat_x, repeat_y);
        this.drawPattern(this.cutOpenPattern, true);
        
        this.statusElement.textContent = `Deployed pattern at ${(openingAngleRad * 180/Math.PI).toFixed(1)}°`;
    }


    drawPatternOnCanvas(pattern, showUnfriendlyVertices, canvas) {
        if (!pattern || !pattern.verts || !pattern.faces) return;
        let ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        
        const verts = pattern.verts;
        const faces = pattern.faces;
        const friendly = pattern.friendly || [];
        
        // Calculate bounding box
        let minX = Infinity, minY = Infinity;
        let maxX = -Infinity, maxY = -Infinity;
        
        verts.forEach(vert => {
            minX = Math.min(minX, vert[0]);
            minY = Math.min(minY, vert[1]);
            maxX = Math.max(maxX, vert[0]);
            maxY = Math.max(maxY, vert[1]);
        });
        
        // Calculate scale and translation
        const width = maxX - minX;
        const height = maxY - minY;
        const scaleX = (canvas.width - 40) / width;
        const scaleY = (canvas.height - 40) / height;
        const scale = Math.min(scaleX, scaleY);
        
        const centerX = (minX + maxX) / 2;
        const centerY = (minY + maxY) / 2;
        const translateX = canvas.width / 2 - centerX * scale;
        const translateY = canvas.height / 2 - centerY * scale;
        
        // Draw faces
        faces.forEach((face, faceIndex) => {
            ctx.beginPath();
            const firstVert = verts[face[0]];
            ctx.moveTo(
                (firstVert[0] * scale + translateX) * this.zoom + this.panX,
                (firstVert[1] * scale + translateY) * this.zoom + this.panY
            );
            
            for (let i = 1; i < face.length; i++) {
                const vert = verts[face[i]];
                ctx.lineTo(
                    (vert[0] * scale + translateX) * this.zoom + this.panX,
                    (vert[1] * scale + translateY) * this.zoom + this.panY
                );
            }
            
            ctx.closePath();
            
            // Fill face
            ctx.fillStyle = `rgba(${faceIndex * 50 % 255}, ${faceIndex * 70 % 255}, ${faceIndex * 90 % 255}, 0.3)`;
            ctx.fill();
            
            // Stroke outline
            ctx.strokeStyle = 'black';
            ctx.lineWidth = 1;
            ctx.stroke();
        });
        
        // Draw unfriendly vertices if requested
        if (showUnfriendlyVertices && friendly.length > 0) {
            verts.forEach((vert, index) => {
                if (friendly[index] === 1) {
                    ctx.beginPath();
                    ctx.arc(
                        (vert[0] * scale + translateX) * this.zoom + this.panX,
                        (vert[1] * scale + translateY) * this.zoom + this.panY,
                        3, 0, Math.PI * 2
                    );
                    ctx.fillStyle = 'rgba(255, 0, 0, 0.8)';
                    ctx.fill();
                    ctx.strokeStyle = 'black';
                    ctx.lineWidth = 1;
                    ctx.stroke();
                }
            });
        }
    }

    drawPattern(pattern, showUnfriendlyVertices) {
      this.drawPatternOnCanvas(pattern, showUnfriendlyVertices, this.canvas);
    }

    handleZoom(deltaY, mouseX, mouseY) {
        const zoomFactor = deltaY > 0 ? 0.9 : 1.1;
        const newZoom = this.zoom * zoomFactor;
        
        // Limit zoom range
        if (newZoom < 0.1 || newZoom > 10) return;
        
        // Calculate mouse position in canvas coordinates
        const mouseCanvasX = (mouseX - this.panX) / this.zoom;
        const mouseCanvasY = (mouseY - this.panY) / this.zoom;
        
        // Update zoom
        this.zoom = newZoom;
        
        // Adjust pan to zoom towards mouse position
        this.panX = mouseX - mouseCanvasX * this.zoom;
        this.panY = mouseY - mouseCanvasY * this.zoom;
        
        this.drawCurrentPattern();
    }

    drawCurrentPattern() {
        if (this.isCutOpen && this.cutOpenPattern) {
            this.drawPattern(this.cutOpenPattern, true);
        } else if (this.basePattern) {
            this.drawPattern(this.cutOpenPattern || this.basePattern, false);
        }
    }

    updatePreview() {
        this.drawCurrentPattern();
    }

    // Public method to check if pattern is valid
    hasValidPattern() {
        return this.basePattern !== null;
    }

    // Add a proper reset method to the PreviewModule class
    reset() {
        this.basePattern = null;
        this.cutOpenPattern = null;
        this.unitPattern = null;
        this.isCutOpen = false;
        
        // Reset zoom and pan state
        this.zoom = 1;
        this.panX = 0;
        this.panY = 0;
        this.autoPlayToggle.checked = true;
        
        // Clear the canvas
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        
        // Reset status message
        this.statusElement.textContent = 'deployment-unfriendly vtx will be highlighted';
        
        // Hide the angle control
        const angleControl = document.querySelector('#problem3Canvas .angle-control');
        if (angleControl) {
            angleControl.style.display = 'none';
        }
        
        // Reset the opening angle slider
        const openingAngleSlider = document.getElementById('openingAngleSlider');
        const angleValue = document.getElementById('angleValue');
        if (openingAngleSlider && angleValue) {
            openingAngleSlider.value = 0;
            angleValue.textContent = '0°';
        }
    }
}