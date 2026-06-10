// js/design-module.js - Parallel line drawing functionality
export class DesignModule {
  constructor(canvasId, statusId) {
    this.canvas = document.getElementById(canvasId);
    this.ctx = this.canvas.getContext('2d');
    this.statusElement = document.getElementById(statusId);

    this.groups = [];
    this.currentGroup = null;
    this.currentGroupIndex = 1;
    this.gridSize = 40;
    this.gridType = 'quad';
    this.draggedPoint = null;
    this.isDragging = false;

    this.initCanvas();
    this.setupEventListeners();
  }

  initCanvas() {
    this.drawGrid();
    this.statusElement.textContent = `Click point A for group ${this.currentGroupIndex}`;
    this.currentGroup = {
      points: [],
      color: this.getGroupColor(this.currentGroupIndex)
    };
  }

  getGroupColor(index) {
    const colors = [
      'rgba(255, 0, 0, 0.8)',    // red
      'rgba(0, 128, 0, 0.8)',    // green
      'rgba(0, 0, 255, 0.8)',     // blue
      'rgba(255, 165, 0, 0.8)',   // orange
      'rgba(128, 0, 128, 0.8)',   // purple
      'rgba(0, 128, 128, 0.8)'    // teal
    ];
    return colors[(index - 1) % colors.length];
  }

  // Grid drawing functions (quad, tri, hex)
  drawGrid() {
    this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

    switch (this.gridType) {
      case 'quad': this.drawQuadGrid(); break;
      case 'triangle': this.drawTriGrid(); break;
      case 'hex': this.drawHexGrid(); break;
    }
  }

  drawQuadGrid() {
    this.ctx.strokeStyle = '#e0e0e0';
    this.ctx.lineWidth = 1;

    // Vertical lines
    for (let x = 0; x <= this.canvas.width; x += this.gridSize) {
      this.ctx.beginPath();
      this.ctx.moveTo(x, 0);
      this.ctx.lineTo(x, this.canvas.height);
      this.ctx.stroke();
    }

    // Horizontal lines
    for (let y = 0; y <= this.canvas.height; y += this.gridSize) {
      this.ctx.beginPath();
      this.ctx.moveTo(0, y);
      this.ctx.lineTo(this.canvas.width, y);
      this.ctx.stroke();
    }
  }

  drawTriGrid() {
    this.ctx.strokeStyle = '#e0e0e0';
    this.ctx.lineWidth = 1;
    const angle = Math.PI / 3;
    const h = this.gridSize * Math.sqrt(3) / 2;

    // Horizontal lines
    for (let y = 0; y <= this.canvas.height; y += h) {
      this.ctx.beginPath();
      this.ctx.moveTo(0, y);
      this.ctx.lineTo(this.canvas.width, y);
      this.ctx.stroke();
    }

    // 60° lines
    let diagonalLength = Math.max(this.canvas.width, this.canvas.height) * 2;
    diagonalLength = Math.ceil(diagonalLength / this.gridSize) * this.gridSize;
    for (let x = -diagonalLength; x <= this.canvas.width + diagonalLength; x += this.gridSize) {
      this.ctx.beginPath();
      this.ctx.moveTo(x, 0);
      this.ctx.lineTo(x - this.canvas.height / Math.tan(angle), this.canvas.height);
      this.ctx.stroke();
    }

    // -60° lines
    for (let x = -diagonalLength; x <= this.canvas.width + diagonalLength; x += this.gridSize) {
      this.ctx.beginPath();
      this.ctx.moveTo(x, 0);
      this.ctx.lineTo(x + this.canvas.height / Math.tan(angle), this.canvas.height);
      this.ctx.stroke();
    }
  }

  drawHexGrid() {
    this.ctx.strokeStyle = '#e0e0e0';
    this.ctx.lineWidth = 1;
    const size = this.gridSize;
    const hexHeight = size * 2;
    const hexWidth = Math.sqrt(3) * size;
    const vertDist = size * 1.5;
    const horizDist = hexWidth;

    const extendedX = this.canvas.width + 2 * hexWidth;
    const extendedY = this.canvas.height + 2 * hexHeight;

    for (let y = -hexHeight; y < extendedY; y += vertDist) {
      for (let x = -hexWidth; x < extendedX; x += horizDist) {
        const xOffset = (Math.floor(y / vertDist) % 2) * (hexWidth / 2);
        this.drawSingleHexagon(x + xOffset, y, size);
      }
    }
  }

  drawSingleHexagon(x, y, size) {
    this.ctx.beginPath();
    for (let i = 0; i < 6; i++) {
      const angle = Math.PI / 3 * i - Math.PI / 6;
      const px = x + size * Math.cos(angle);
      const py = y + size * Math.sin(angle);
      if (i === 0) {
        this.ctx.moveTo(px, py);
      } else {
        this.ctx.lineTo(px, py);
      }
    }
    this.ctx.closePath();
    this.ctx.stroke();
  }

  // Snapping functions
  snapToGrid(x, y) {
    switch (this.gridType) {
      case 'quad': return this.snapQuad(x, y);
      case 'triangle': return this.snapTri(x, y);
      case 'hex': return this.snapHex(x, y);
      default: return { x, y };
    }
  }

  snapQuad(x, y) {
    return {
      x: Math.round(x / this.gridSize) * this.gridSize,
      y: Math.round(y / this.gridSize) * this.gridSize
    };
  }

  snapTri(x, y) {
    const h = this.gridSize * Math.sqrt(3) / 2;
    const row = Math.round(y / h);
    const ySnapped = row * h;
    const xOffset = (row % 2) * this.gridSize / 2;
    const col = Math.round((x - xOffset) / this.gridSize);
    const xSnapped = col * this.gridSize + xOffset;

    return { x: xSnapped, y: ySnapped };
  }

  snapHex(x, y) {
    const size = this.gridSize;
    const q = (x * Math.sqrt(3) / 3 - y / 3) / size;
    const r = y * 2 / 3 / size;

    let rq = Math.round(q);
    let rr = Math.round(r);
    const rs = Math.round(-q - r);

    const q_diff = Math.abs(rq - q);
    const r_diff = Math.abs(rr - r);
    const s_diff = Math.abs(rs - (-q - r));

    if (q_diff > r_diff && q_diff > s_diff) {
      rq = -rr - rs;
    } else if (r_diff > s_diff) {
      rr = -rq - rs;
    }

    return {
      x: size * Math.sqrt(3) * (rq + rr / 2),
      y: size * 3 / 2 * rr
    };
  }

  // Drawing functions
  drawPoints() {
    // Draw all groups' points
    this.groups.forEach((group, groupIdx) => {
      group.points.forEach((point, pointIdx) => {
        this.ctx.beginPath();
        this.ctx.arc(point.x, point.y, 6, 0, Math.PI * 2);
        this.ctx.fillStyle = group.color;
        this.ctx.fill();
        this.ctx.strokeStyle = 'black';
        this.ctx.stroke();

        // Label points
        this.ctx.fillStyle = 'black';
        this.ctx.font = '14px Arial';
        const label = ['A', 'B', 'C'][pointIdx] + (groupIdx + 1);
        this.ctx.fillText(label, point.x + 10, point.y - 10);
      });
    });

    // Draw current group points
    if (this.currentGroup && this.currentGroup.points.length > 0) {
      this.currentGroup.points.forEach((point, pointIdx) => {
        this.ctx.beginPath();
        this.ctx.arc(point.x, point.y, 6, 0, Math.PI * 2);
        this.ctx.fillStyle = this.currentGroup.color;
        this.ctx.fill();
        this.ctx.strokeStyle = 'black';
        this.ctx.stroke();

        // Label points
        this.ctx.fillStyle = 'black';
        this.ctx.font = '14px Arial';
        const label = ['A', 'B', 'C'][pointIdx] + this.currentGroupIndex;
        this.ctx.fillText(label, point.x + 10, point.y - 10);
      });
    }
  }

  drawParallelLines() {
    // Draw completed groups
    this.groups.forEach(group => {
      if (group.points.length < 3) return;
      this.drawSingleGroupLines(group);
    });

    // Draw current group
    if (this.currentGroup && this.currentGroup.points.length >= 2) {
      this.drawSingleGroupLines(this.currentGroup, true);
    }
  }

  drawSingleGroupLines(group, isCurrent = false) {
    const A = group.points[0];
    const B = group.points[1];

    // If we don't have at least 2 points, return
    if (!A || !B) return;

    const C = group.points[2] || group.points[1];

    // Calculate line direction vector
    const dir = { x: B.x - A.x, y: B.y - A.y };

    // Handle case where A and B are the same point
    if (Math.abs(dir.x) < 0.001 && Math.abs(dir.y) < 0.001) return;

    // Calculate normal vector
    const normal = { x: -dir.y, y: dir.x };
    const normalLength = Math.sqrt(normal.x * normal.x + normal.y * normal.y);

    // Handle case where normal length is zero
    if (normalLength < 0.001) return;

    const unitNormal = { x: normal.x / normalLength, y: normal.y / normalLength };

    // Calculate offset distance
    const offsetDistance = this.distanceToLine(A, B, C);
    const sign = ((C.x - A.x) * normal.x + (C.y - A.y) * normal.y) > 0 ? 1 : -1;
    const signedOffset = offsetDistance * sign;

    // Draw the main line AB
    this.ctx.strokeStyle = 'black';
    this.ctx.lineWidth = isCurrent ? 1 : 2;
    this.ctx.beginPath();

    // Extend line to canvas edges
    if (Math.abs(dir.x) < 0.001) {
      // Vertical line
      this.ctx.moveTo(A.x, 0);
      this.ctx.lineTo(A.x, this.canvas.height);
    } else if (Math.abs(dir.y) < 0.001) {
      // Horizontal line
      this.ctx.moveTo(0, A.y);
      this.ctx.lineTo(this.canvas.width, A.y);
    } else {
      // Diagonal line - extend to canvas boundaries
      const slope = dir.y / dir.x;
      const yAtX0 = A.y - slope * A.x;
      const yAtXMax = A.y + slope * (this.canvas.width - A.x);

      this.ctx.moveTo(0, yAtX0);
      this.ctx.lineTo(this.canvas.width, yAtXMax);
    }
    this.ctx.stroke();

    // Draw parallel lines if we have all 3 points
    if (group.points.length >= 3 && offsetDistance > 0.001) {
      this.ctx.strokeStyle = group.color.replace('0.8', '0.4');
      this.ctx.lineWidth = 1;

      const canvasDiagonal = Math.sqrt(
        this.canvas.width * this.canvas.width +
        this.canvas.height * this.canvas.height
      );
      const numLines = Math.ceil(canvasDiagonal / offsetDistance) + 2;

      for (let i = -numLines; i <= numLines; i++) {
        if (i === 0) continue;

        const offset = i * signedOffset;
        const offsetX = unitNormal.x * offset;
        const offsetY = unitNormal.y * offset;

        this.ctx.beginPath();

        if (Math.abs(dir.x) < 0.001) {
          // Vertical line
          this.ctx.moveTo(A.x + offsetX, 0);
          this.ctx.lineTo(A.x + offsetX, this.canvas.height);
        } else if (Math.abs(dir.y) < 0.001) {
          // Horizontal line
          this.ctx.moveTo(0, A.y + offsetY);
          this.ctx.lineTo(this.canvas.width, A.y + offsetY);
        } else {
          // Diagonal line
          const slope = dir.y / dir.x;
          const yAtX0 = (A.y + offsetY) - slope * (A.x + offsetX);
          const yAtXMax = yAtX0 + slope * this.canvas.width;

          this.ctx.moveTo(0, yAtX0);
          this.ctx.lineTo(this.canvas.width, yAtXMax);
        }
        this.ctx.stroke();
      }
    }
  }

  distanceToLine(A, B, C) {
    // Calculate the distance from point C to line AB
    const numerator = Math.abs(
      (B.y - A.y) * C.x -
      (B.x - A.x) * C.y +
      B.x * A.y -
      B.y * A.x
    );
    const denominator = Math.sqrt(
      Math.pow(B.y - A.y, 2) +
      Math.pow(B.x - A.x, 2)
    );
    return numerator / denominator;
  }

  findPointAt(x, y) {
    // Check current group
    if (this.currentGroup && this.currentGroup.points.length > 0) {
      for (let p = 0; p < this.currentGroup.points.length; p++) {
        const point = this.currentGroup.points[p];
        const distance = Math.sqrt(Math.pow(x - point.x, 2) + Math.pow(y - point.y, 2));
        if (distance <= 10) {
          return { groupIndex: -1, pointIndex: p, point: point };
        }
      }
    }

    // Check completed groups
    for (let g = 0; g < this.groups.length; g++) {
      for (let p = 0; p < this.groups[g].points.length; p++) {
        const point = this.groups[g].points[p];
        const distance = Math.sqrt(Math.pow(x - point.x, 2) + Math.pow(y - point.y, 2));
        if (distance <= 10) {
          return { groupIndex: g, pointIndex: p, point: point };
        }
      }
    }
    return null;
  }

  redrawCanvas() {
    this.drawGrid();
    this.drawParallelLines();
    this.drawPoints();
  }

  setupEventListeners() {
    // Click to add points
    this.canvas.addEventListener('click', (e) => {
      if (this.isDragging) {
        this.isDragging = false;
        return;
      }

      if (!this.currentGroup) {
        this.currentGroup = {
          points: [],
          color: this.getGroupColor(this.currentGroupIndex)
        };
      }

      if (this.currentGroup.points.length < 3) {
        const rect = this.canvas.getBoundingClientRect();
        const rawX = e.clientX - rect.left;
        const rawY = e.clientY - rect.top;
        const { x, y } = this.snapToGrid(rawX, rawY);

        this.currentGroup.points.push({ x, y });

        if (this.currentGroup.points.length === 1) {
          this.statusElement.textContent = `Click point B for group ${this.currentGroupIndex}`;
        } else if (this.currentGroup.points.length === 2) {
          this.statusElement.textContent = `Click point C for group ${this.currentGroupIndex}`;
        } else if (this.currentGroup.points.length === 3) {
          this.groups.push(this.currentGroup);
          this.currentGroup = null;
          this.currentGroupIndex++;
          this.statusElement.textContent = `Group ${this.groups.length} complete. Continue clicking to add another group or drag points to adjust.`;
        }

        this.redrawCanvas();
      }
    });

    // Mouse down for dragging
    this.canvas.addEventListener('mousedown', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;

      const foundPoint = this.findPointAt(x, y);
      if (foundPoint) {
        this.draggedPoint = foundPoint;
        this.isDragging = true;
        this.canvas.style.cursor = 'grabbing';
      }
    });

    // Mouse move for dragging
    this.canvas.addEventListener('mousemove', (e) => {
      if (this.draggedPoint) {
        const rect = this.canvas.getBoundingClientRect();
        const rawX = e.clientX - rect.left;
        const rawY = e.clientY - rect.top;
        const { x, y } = this.snapToGrid(rawX, rawY);

        this.draggedPoint.point.x = x;
        this.draggedPoint.point.y = y;
        this.redrawCanvas();
      } else {
        // Change cursor when hovering over a point
        const rect = this.canvas.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        const foundPoint = this.findPointAt(x, y);
        this.canvas.style.cursor = foundPoint ? 'grab' : 'crosshair';
      }
    });

    // Mouse up after dragging
    this.canvas.addEventListener('mouseup', () => {
      if (this.draggedPoint) {
        this.draggedPoint = null;
        this.canvas.style.cursor = 'crosshair';
      }
    });

    // Mouse leave while dragging
    this.canvas.addEventListener('mouseleave', () => {
      if (this.draggedPoint) {
        this.draggedPoint = null;
        this.canvas.style.cursor = 'crosshair';
      }
    });

    // document.getElementById('savePatternBtn').addEventListener('click', () => {
    //   const data = JSON.stringify(this.groups);
    //   const blob = new Blob([data], { type: 'application/json' });
    //   const url = URL.createObjectURL(blob);
    //   const a = document.createElement('a');
    //   a.href = url;
    //   const exportName = 'parallel_lines_groups';
    //   a.download = `${exportName}.json`;
    //   document.body.appendChild(a);
    //   a.click();
    //   document.body.removeChild(a);
    //   URL.revokeObjectURL(url);
    // });

    document.getElementById('loadPatternBtn').addEventListener('click', () => {
      const fileInput = document.createElement('input');
      fileInput.type = 'file';
      fileInput.accept = '.json';
      fileInput.onchange = (e) => {
        const file = e.target.files[0];
        const reader = new FileReader();
        reader.onload = (e) => {
          const data = JSON.parse(e.target.result);
          this.reset();
          this.groups = data['groups'];
          this.currentGroup = null;
          this.currentGroupIndex = this.groups.length + 1;
          this.redrawCanvas();
          this.triggerPreviewUpdate();
          if (window.previewModule) {
            previewModule.generateBasePattern();
            previewModule.updatePreview();
          }
        };
        reader.readAsText(file);
      };
      fileInput.click();
    });

    document.getElementById('patternPresetSelect').addEventListener('change', async (e) => {
      this.loadPresetPattern(e.target.value);
    });
  }

  async loadPresetPatternGroups(preset) {
    if (!preset)
      return [];
    try {
      const presetJson = await fetch('data/patterns/' + preset + '.json');
      const importedGroups = await presetJson.json();
      return importedGroups['groups'];
    } catch (error) {
      console.error('Error loading preset: ' + error);
      return [];
    }
  }
  async loadPresetPattern(preset) {
    if (!preset)
      return;
    try {
      const presetJson = await fetch('data/patterns/' + preset + '.json');
      const data = await presetJson.json();
      this.reset();
      this.groups = data['groups'];
      this.gridSize = data['gridSize'] || 40;
      this.gridType = data['gridType'] || 'quad';
      document.getElementById('gridType').value = this.gridType;
      document.getElementById('gridSize').value = this.gridSize;
      document.getElementById('gridSizeValue').textContent = this.gridSize;
      this.drawGrid();
      this.currentGroup = null;
      this.currentGroupIndex = this.groups.length + 1;
      this.redrawCanvas();
      this.triggerPreviewUpdate();
    } catch (error) {
      console.error('Error loading preset: ' + error);
    }
  }

  // Public methods
  reset() {
    this.groups = [];
    this.currentGroup = null;
    this.currentGroupIndex = 1;
    this.initCanvas();
    this.redrawCanvas();
    this.triggerPreviewUpdate();
  }

  enable() {
    this.canvas.style.pointerEvents = 'auto';
    this.canvas.style.opacity = '1';
    this.canvas.style.cursor = 'crosshair';
  }

  disable() {
    this.canvas.style.pointerEvents = 'none';
    this.canvas.style.opacity = '0.5';
    this.canvas.style.cursor = 'default';
  }

  findIntersection(line1, line2) {
    const x1 = line1.x1, y1 = line1.y1;
    const x2 = line1.x2, y2 = line1.y2;
    const x3 = line2.x1, y3 = line2.y1;
    const x4 = line2.x2, y4 = line2.y2;

    const denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (denominator === 0) return null; // Lines are parallel

    const t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denominator;
    const u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denominator;

    // Check if intersection is within both line segments
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
      return {
        x: x1 + t * (x2 - x1),
        y: y1 + t * (y2 - y1)
      };
    }
    return null;
  }

  isPointInBounds(x, y) {
    const margin = 100; // Allow points slightly outside the canvas
    return x >= -margin && x <= this.canvas.width + margin &&
      y >= -margin && y <= this.canvas.height + margin;
  }

  generateGraph() {
    const vertices = [];
    const edges = [];
    const vertexMap = new Map();
    const lineIntersections = new Map();

    // Collect all lines from all groups
    const allLines = [];
    this.groups.forEach(group => {
      if (group.points.length < 3) return;

      const A = group.points[0];
      const B = group.points[1];
      const C = group.points[2];

      // Calculate line direction and normal
      const dir = { x: B.x - A.x, y: B.y - A.y };
      const normal = { x: -dir.y, y: dir.x };
      const normalLength = Math.sqrt(normal.x * normal.x + normal.y * normal.y);

      if (normalLength < 0.001) return; // Skip if normal is zero

      const unitNormal = {
        x: normal.x / normalLength,
        y: normal.y / normalLength
      };

      const offsetDistance = this.distanceToLine(A, B, C);
      const sign = ((C.x - A.x) * normal.x + (C.y - A.y) * normal.y) > 0 ? 1 : -1;
      const signedOffset = offsetDistance * sign;

      // Add main line
      if (Math.abs(dir.x) < 0.001) {
        // Vertical line
        allLines.push({
          x1: A.x, y1: 0,
          x2: A.x, y2: this.canvas.height,
          group: group
        });
      } else if (Math.abs(dir.y) < 0.001) {
        // Horizontal line
        allLines.push({
          x1: 0, y1: A.y,
          x2: this.canvas.width, y2: A.y,
          group: group
        });
      } else {
        // Diagonal line
        const slope = dir.y / dir.x;
        const yAtX0 = A.y - slope * A.x;
        const yAtXMax = A.y + slope * (this.canvas.width - A.x);

        allLines.push({
          x1: 0, y1: yAtX0,
          x2: this.canvas.width, y2: yAtXMax,
          group: group
        });
      }

      // Add parallel lines if we have all 3 points
      if (group.points.length >= 3 && offsetDistance > 0.001) {
        const canvasDiagonal = Math.sqrt(
          this.canvas.width * this.canvas.width +
          this.canvas.height * this.canvas.height
        );
        const numLines = Math.ceil(canvasDiagonal / offsetDistance) + 2;

        for (let i = -numLines; i <= numLines; i++) {
          if (i === 0) continue;

          const offset = i * signedOffset;
          const offsetX = unitNormal.x * offset;
          const offsetY = unitNormal.y * offset;

          if (Math.abs(dir.x) < 0.001) {
            // Vertical parallel lines
            allLines.push({
              x1: A.x + offsetX, y1: 0,
              x2: A.x + offsetX, y2: this.canvas.height,
              group: group
            });
          } else if (Math.abs(dir.y) < 0.001) {
            // Horizontal parallel lines
            allLines.push({
              x1: 0, y1: A.y + offsetY,
              x2: this.canvas.width, y2: A.y + offsetY,
              group: group
            });
          } else {
            // Diagonal parallel lines
            const slope = dir.y / dir.x;
            const yAtX0 = (A.y + offsetY) - slope * (A.x + offsetX);
            const yAtXMax = yAtX0 + slope * this.canvas.width;

            allLines.push({
              x1: 0, y1: yAtX0,
              x2: this.canvas.width, y2: yAtXMax,
              group: group
            });
          }
        }
      }
    });

    // Initialize lineIntersections map
    allLines.forEach((_, index) => {
      lineIntersections.set(index, []);
    });

    // Find all intersections
    for (let i = 0; i < allLines.length; i++) {
      for (let j = i + 1; j < allLines.length; j++) {
        const intersection = this.findIntersection(allLines[i], allLines[j]);
        if (intersection && this.isPointInBounds(intersection.x, intersection.y)) {
          lineIntersections.get(i).push(intersection);
          lineIntersections.get(j).push(intersection);
        }
      }
    }

    // Sort intersections along each line and create edges
    allLines.forEach((line, lineIndex) => {
      const intersections = lineIntersections.get(lineIndex);

      if (intersections.length === 0) return;

      // Sort intersections along the line
      intersections.sort((a, b) => {
        // For vertical lines, sort by y-coordinate
        if (line.x1 === line.x2) {
          return a.y - b.y;
        }
        // For non-vertical lines, sort by x-coordinate
        return a.x - b.x;
      });

      // Create edges between consecutive points
      for (let i = 0; i < intersections.length - 1; i++) {
        const p1 = intersections[i];
        const p2 = intersections[i + 1];

        // Add vertices if they don't exist
        const key1 = `${p1.x},${p1.y}`;
        const key2 = `${p2.x},${p2.y}`;

        if (!vertexMap.has(key1)) {
          vertexMap.set(key1, vertices.length);
          vertices.push(p1);
        }
        if (!vertexMap.has(key2)) {
          vertexMap.set(key2, vertices.length);
          vertices.push(p2);
        }

        // Add edge between consecutive points
        edges.push({
          v1: vertexMap.get(key1),
          v2: vertexMap.get(key2),
          group: line.group
        });
      }
    });

    return { vertices, edges };
  }

  // Add this method to check deployment-friendly vertices
  checkDeploymentFriendlyVertices() {
    const graph = this.generateGraph();
    const unfriendlyVertices = [];

    // For each vertex, check if it's deployment-friendly
    // (sum of directed edges ending at the vertex should be zero)
    graph.vertices.forEach((vertex, vertexIndex) => {
      // Find all edges that end at this vertex
      const incomingEdges = graph.edges.filter(edge => {
        const edgeEnd = graph.vertices[edge.v2];
        return Math.abs(edgeEnd.x - vertex.x) < 0.001 &&
          Math.abs(edgeEnd.y - vertex.y) < 0.001;
      });

      if (incomingEdges.length > 0) {
        // Calculate vector sum of incoming edges
        let sumX = 0;
        let sumY = 0;

        incomingEdges.forEach(edge => {
          const start = graph.vertices[edge.v1];
          const end = graph.vertices[edge.v2];
          sumX += (end.x - start.x);
          sumY += (end.y - start.y);
        });

        // Check if sum is not zero (with some tolerance)
        if (Math.abs(sumX) > 0.001 || Math.abs(sumY) > 0.001) {
          unfriendlyVertices.push(vertexIndex);
        }
      }
    });

    return unfriendlyVertices;
  }

  // Add trigger method for preview updates
  triggerPreviewUpdate() {
    if (window.previewModule) {
      previewModule.updatePreview();
    }
  }

  // Modify redrawCanvas to trigger preview updates
  redrawCanvas() {
    this.drawGrid();
    this.drawParallelLines();
    this.drawPoints();
    this.triggerPreviewUpdate(); // Add this line
  }

  // Also add triggerPreviewUpdate() to other methods that change the design
  reset() {
    this.groups = [];
    this.currentGroup = null;
    this.currentGroupIndex = 1;
    this.initCanvas();
    this.triggerPreviewUpdate(); // Add this
  }

  setGridType(type) {
    if (this.gridType === type) return;
    this.reset();
    this.gridType = type;
    this.redrawCanvas();
    this.triggerPreviewUpdate(); // Add this
  }

  setGridSize(size) {
    console.log("Setting grid size to:", size);
    const previousSize = this.gridSize || 1;
    this.gridSize = size;

    // Scale all existing points by the ratio of new/old grid size,
    // then snap once to the new grid to keep them aligned.
    const scale = size / previousSize;
    this.groups.forEach(group => {
      group.points = group.points.map(p => {
        const scaledX = p.x * scale;
        const scaledY = p.y * scale;
        const snapped = this.snapToGrid(scaledX, scaledY);
        return { x: snapped.x, y: snapped.y };
      });
    });
    if (this.currentGroup && this.currentGroup.points) {
      this.currentGroup.points = this.currentGroup.points.map(p => {
        const scaledX = p.x * scale;
        const scaledY = p.y * scale;
        const snapped = this.snapToGrid(scaledX, scaledY);
        return { x: snapped.x, y: snapped.y };
      });
    }

    this.redrawCanvas();
    this.triggerPreviewUpdate(); // Add this
  }

  generateDeploymentVectors() {
    const graph = this.generateGraph();
    const deploymentVectors = [];

    // For each edge, calculate the deployment direction
    graph.edges.forEach(edge => {
      const v1 = graph.vertices[edge.v1];
      const v2 = graph.vertices[edge.v2];

      // Calculate edge center and normal direction for deployment
      const center = {
        x: (v1.x + v2.x) / 2,
        y: (v1.y + v2.y) / 2
      };

      // Calculate normal vector (perpendicular to edge)
      const dx = v2.x - v1.x;
      const dy = v2.y - v1.y;
      const length = Math.sqrt(dx * dx + dy * dy);

      if (length > 0.001) {
        const normal = {
          x: -dy / length,
          y: dx / length
        };

        deploymentVectors.push({
          center: center,
          normal: normal,
          group: edge.group
        });
      }
    });

    return deploymentVectors;
  }

}