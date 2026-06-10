import { DesignModule } from './design-module.js';
import { PreviewModule } from './preview-module.js';
import { AlgorithmModule } from './algorithm-module.js';

document.addEventListener('DOMContentLoaded', function() {
    
    // Set canvas dimensions for design module
    const designCanvas = document.getElementById('designCanvas');
    const canvasSize = Math.min(designCanvas.parentElement.clientWidth, 400);
    designCanvas.width = canvasSize;
    designCanvas.height = canvasSize;
    
    // Create the design and preview modules
    window.designModule = new DesignModule('designCanvas', 'designStatus');
    window.previewModule = new PreviewModule('previewCanvas', 'previewStatus', window.designModule);
    
    // Initially disable the design module (it's inside a closed dropdown)
    designModule.disable();
    designModule.redrawCanvas();

    // --- Problem 3: Design Canvas Event Listeners ---
    const resetBtn = document.getElementById('resetBtn');
    const gridSizeSlider = document.getElementById('gridSize');
    const gridSizeValue = document.getElementById('gridSizeValue');
    const gridTypeSelect = document.getElementById('gridType');
    const generatePatternBtn = document.getElementById('generatePatternBtn');
    const repeatXInput = document.getElementById('repeatX');
    const repeatYInput = document.getElementById('repeatY');
    const savePatternBtn = document.getElementById('savePatternBtn')

    if (resetBtn) {
        resetBtn.addEventListener('click', () => {
            // Reset both modules
            designModule.reset();
            previewModule.reset();
            
            // You might also need to reset any other state variables
            // For example, reset the grid if needed
            const gridSize = document.getElementById('gridSize');
            const gridType = document.getElementById('gridType');
            if (gridSize && gridType) {
                gridSize.value = 40;
                gridType.value = 'quad';
            }
            
            // Reset repeat values
            const repeatX = document.getElementById('repeatX');
            const repeatY = document.getElementById('repeatY');
            if (repeatX && repeatY) {
                repeatX.value = 3;
                repeatY.value = 3;
            }
            
        });
    }

    if (gridSizeSlider) {
        gridSizeSlider.addEventListener('input', function(e) {
            const size = parseInt(e.target.value);
            gridSizeValue.textContent = size;
            if (window.designModule) designModule.setGridSize(size);
        });
    }

    if (gridTypeSelect) {
        gridTypeSelect.addEventListener('change', function(e) {
            if (window.designModule) designModule.setGridType(e.target.value);
        });
    }

    if (generatePatternBtn) {
        generatePatternBtn.addEventListener('click', function() {
            if (window.previewModule) {
                previewModule.generateBasePattern();
                previewModule.updatePreview();
            }
        });
    }

    if (repeatXInput) {
        repeatXInput.addEventListener('input', () => {
            if (window.previewModule && previewModule.isCutOpen && previewModule.basePattern) {
                previewModule.cut_open_pattern();
            }
        });
    }

    if (repeatYInput) {
        repeatYInput.addEventListener('input', () => {
            if (window.previewModule && previewModule.isCutOpen && previewModule.basePattern) {
                previewModule.cut_open_pattern();
            }
        });
    }


    // Add event listener for the save button
    if (savePatternBtn) {
        savePatternBtn.addEventListener('click', savePatternWithPrompt);
    }

    // --- Problem 4: Algorithm Module Listeners ---

    // Initialize Algorithm Module with Three.js container
    window.algorithmModule = new AlgorithmModule();
    
    // Initialize Three.js with the meshViewer container
    window.algorithmModule.init();
    // Set up mesh loading event listeners
    const loadMeshBtn = document.getElementById('loadMeshBtn');
    const meshFileInput = document.getElementById('meshFileInput');
    const toggleWireframeBtn = document.getElementById('toggleWireframeBtn');

    if (loadMeshBtn && meshFileInput) {
        loadMeshBtn.addEventListener('click', () => {
            meshFileInput.click();
        });
        
        meshFileInput.addEventListener('change', (e) => {
            window.algorithmModule.handleMeshFileSelect(e);
        });
    }
    
    if (toggleWireframeBtn) {
        toggleWireframeBtn.addEventListener('click', () => {
            window.algorithmModule.toggleWireframe();
        });
    }

    // --- Dropdown Toggle Logic ---
    // Enable/disable design canvas based on which problem dropdown is open
    document.querySelectorAll('.solution-toggle-btn').forEach(button => {
        button.addEventListener('click', function(e) {
            e.stopPropagation();
            const canvas = this.closest('.problem-container').querySelector('.dropdown-canvas');
            const isOpening = !canvas.classList.contains('active');
            const isProblem3 = this.closest('#problem3');

            // Close all other dropdowns first
            document.querySelectorAll('.dropdown-canvas.active').forEach(otherCanvas => {
                if (otherCanvas !== canvas) {
                    otherCanvas.classList.remove('active');
                    const otherButton = otherCanvas.closest('.problem-container').querySelector('.solution-toggle-btn');
                    otherButton.textContent = 'show our solution';
                    // Disable design module if another dropdown was open
                    if (window.designModule && !isProblem3) {
                        designModule.disable();
                    }
                }
            });

            // Toggle the current dropdown
            canvas.classList.toggle('active');
            this.textContent = canvas.classList.contains('active') ? 'hide our solution' : 'show our solution';

            // Enable design module ONLY if Problem 3 is being opened
            if (window.designModule) {
                if (isProblem3 && isOpening) {
                    // Small timeout to ensure the dropdown is open before enabling
                    setTimeout(() => {
                        designModule.enable();
                        designModule.redrawCanvas();
                    }, 50);
                } else if (isProblem3 && !isOpening) {
                    // Optional: disable when Problem 3 is closed
                    // designModule.disable();
                }
            }
        });
    });

    // --- Algorithm Overview Toggle (Problem 4) ---
    const overviewToggleBtn = document.getElementById('overviewToggleBtn');
    const algorithmOverview = document.getElementById('algorithmOverview');
    
    if (overviewToggleBtn && algorithmOverview) {
        algorithmOverview.classList.add('expanded'); // Start expanded
        overviewToggleBtn.addEventListener('click', function() {
            const isExpanded = algorithmOverview.classList.contains('expanded');
            if (isExpanded) {
                algorithmOverview.classList.remove('expanded');
                algorithmOverview.classList.add('collapsed');
                overviewToggleBtn.textContent = 'show more';
            } else {
                algorithmOverview.classList.remove('collapsed');
                algorithmOverview.classList.add('expanded');
                overviewToggleBtn.textContent = 'show less';
            }
        });
    }

    // --- Window Resize Handler for Design Canvas ---
    window.addEventListener('resize', function() {
        if (window.designModule) {
            const designCanvas = document.getElementById('designCanvas');
            const canvasSize = Math.min(designCanvas.parentElement.clientWidth, 400);
            designCanvas.width = canvasSize;
            designCanvas.height = canvasSize;
            designModule.redrawCanvas();
        }
        
        // Handle Three.js resize
        if (window.algorithmModule) {
            window.algorithmModule.onWindowResize();
        }
    });



    // Optimization weight controls
    const planarityWeight = document.getElementById('planarityWeight');
    const shapeWeight = document.getElementById('shapeWeight');
    const rigidityWeight = document.getElementById('rigidityWeight');
    const fairnessWeight = document.getElementById('fairnessWeight');
    const runOptimizationBtn = document.getElementById('runOptimizationBtn');

    // Store current weights
    let optimizationWeights = {
        planarity: 5,
        shape: 5,
        rigidity: 5,
        fairness: 5
    };

    // Add event listeners for weight inputs
    if (planarityWeight) {
        planarityWeight.addEventListener('input', function(e) {
            optimizationWeights.planarity = parseFloat(e.target.value);
            updateOptimizationStatus();
        });
    }

    if (shapeWeight) {
        shapeWeight.addEventListener('input', function(e) {
            optimizationWeights.shape = parseFloat(e.target.value);
            updateOptimizationStatus();
        });
    }

    if (rigidityWeight) {
        rigidityWeight.addEventListener('input', function(e) {
            optimizationWeights.rigidity = parseFloat(e.target.value);
            updateOptimizationStatus();
        });
    }

    if (fairnessWeight) {
        fairnessWeight.addEventListener('input', function(e) {
            optimizationWeights.fairness = parseFloat(e.target.value);
            updateOptimizationStatus();
        });
    }
});



// Add save pattern functionality
// Replace the savePatternWithPrompt function in app.js with this:
function savePatternWithPrompt() {
    // Check if we have a pattern to save
    if (!window.designModule || window.designModule.groups.length === 0) {
        alert('No pattern to save. Please create a pattern first.');
        return;
    }
    
    // Get filename from user
    const defaultName = 'kirigami_pattern_' + new Date().toISOString().slice(0, 10);
    const filename = prompt('Enter a name for your pattern:', defaultName);
    
    if (filename === null) {
        // User cancelled
        return;
    }
    
    if (filename.trim() === '') {
        alert('Please enter a valid filename.');
        return;
    }
    
    // Get the pattern data (the groups array, grid size, and grid type)
    const patternData = {
        groups: window.designModule.groups,
        gridSize: window.designModule.gridSize,
        gridType: window.designModule.gridType
    };

    // Create a Blob and download link
    const blob = new Blob([JSON.stringify(patternData, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    
    // Create a temporary download link
    const a = document.createElement('a');
    a.href = url;
    a.download = filename.trim() + '.json';
    document.body.appendChild(a);
    a.click();
    
    // Clean up
    setTimeout(() => {
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }, 100);
}

