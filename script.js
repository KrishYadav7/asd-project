function runAnalysis() {
    // 1. Fetch Inputs
    const E = parseFloat(document.getElementById('E').value);
    const rho = parseFloat(document.getElementById('rho').value);
    const A = parseFloat(document.getElementById('A').value);
    const L = parseFloat(document.getElementById('L').value);
    const N = parseInt(document.getElementById('N').value);
    const plotMode = parseInt(document.getElementById('mode').value);

    let logText = `--- STARTING CALCULATION ---\n`;
    logText += `Parameters: E=${E} Pa, rho=${rho} kg/m³, A=${A} m², L=${L} m, Elements N=${N}\n\n`;

    // 2. Element Properties
    const Le = L / N;
    const n_dof = N - 1; // Active degrees of freedom (Nodes 1 to N-1)
    
    const ke = (E * A) / Le;
    const me = (rho * A * Le) / 6;

    logText += `Element Length (Le): ${Le} m\n`;
    logText += `Local Stiffness Multiplier (ke): ${ke.toExponential(4)}\n`;
    logText += `Local Mass Multiplier (me): ${me.toExponential(4)}\n\n`;

    if (N < 2) {
        alert("Number of elements must be at least 2 for Fixed-Fixed.");
        return;
    }

    // 3. Initialize Global Matrices
    let K = numeric.rep([n_dof, n_dof], 0);
    let M = numeric.rep([n_dof, n_dof], 0);

    // 4. Assemble Global Matrices and Apply Boundary Conditions
    // Nodes are 0 to N. Fixed at 0 and N.
    // Loop through elements e = 0 to N-1
    for (let e = 0; e < N; e++) {
        let node1 = e;
        let node2 = e + 1;

        // Map global nodes to reduced active matrix indices (node 1 becomes index 0)
        let i = node1 - 1; 
        let j = node2 - 1; 

        if (i >= 0 && i < n_dof) {
            K[i][i] += ke * 1;
            M[i][i] += me * 2;
        }
        if (j >= 0 && j < n_dof) {
            K[j][j] += ke * 1;
            M[j][j] += me * 2;
        }
        if (i >= 0 && i < n_dof && j >= 0 && j < n_dof) {
            K[i][j] += ke * -1;
            K[j][i] += ke * -1;
            M[i][j] += me * 1;
            M[j][i] += me * 1;
        }
    }

    logText += `--- GLOBAL MATRICES ASSEMBLED ---\n`;
    logText += `Active Degrees of Freedom: ${n_dof}x${n_dof}\n`;
    logText += `(Matrix boundaries applied: First and last nodes removed)\n\n`;

    // 5. Solve Generalized Eigenvalue Problem: (K - w^2 * M) * phi = 0
    // Equivalent to standard problem: (M^-1 * K) * phi = w^2 * phi
    logText += `Inverting Global Mass Matrix [M]...\n`;
    let Minv;
    try {
        Minv = numeric.inv(M);
    } catch(err) {
        document.getElementById('calcLog').innerText = "Matrix inversion failed. Check inputs.";
        return;
    }
    
    logText += `Calculating System Matrix [D] = [M]^-1 * [K]...\n`;
    let D = numeric.dot(Minv, K);
    
    logText += `Solving Eigenvalues (λ = ω²)...\n\n`;
    let eig = numeric.eig(D);
    
    // 6. Extract, Calculate, and Sort Frequencies
    let eigenvalues = eig.lambda.x; // Real parts of eigenvalues
    let eigenvectors = eig.E.x;     // Matrix where columns are eigenvectors

    let frequenciesFEM = [];
    for (let i = 0; i < eigenvalues.length; i++) {
        let lambda = Math.max(0, eigenvalues[i]); // Prevent tiny negative floats from numerical error
        frequenciesFEM.push({
            omega: Math.sqrt(lambda),
            index: i
        });
    }
    
    // Sort ascending by frequency
    frequenciesFEM.sort((a, b) => a.omega - b.omega);

    // 7. Calculate Analytical Frequencies & Generate Table
    const c = Math.sqrt(E / rho);
    let resultsHTML = `
        <table>
            <tr>
                <th>Mode (n)</th>
                <th>Analytical &omega; (rad/s)</th>
                <th>FEM &omega; (rad/s)</th>
                <th>Error (%)</th>
            </tr>`;

    let modesToShow = Math.min(5, n_dof); // Show up to 5 modes in table
    for (let n = 1; n <= modesToShow; n++) {
        let w_ana = (n * Math.PI / L) * c;
        let w_fem = frequenciesFEM[n-1].omega;
        let error = Math.abs((w_fem - w_ana) / w_ana) * 100;

        resultsHTML += `
            <tr>
                <td>${n}</td>
                <td>${w_ana.toFixed(2)}</td>
                <td>${w_fem.toFixed(2)}</td>
                <td>${error.toFixed(4)}%</td>
            </tr>`;
            
        logText += `Mode ${n}: w_ana = ${w_ana.toFixed(2)}, w_fem = ${w_fem.toFixed(2)}, Error = ${error.toFixed(4)}%\n`;
    }
    resultsHTML += `</table>`;
    document.getElementById('resultsTableContainer').innerHTML = resultsHTML;

    // 8. Extract Selected Mode Shape for Plotting
    let modeIndex = plotMode;
    if(modeIndex > n_dof) {
        alert(`Requested mode ${modeIndex} exceeds active DOFs (${n_dof}). Defaulting to Mode 1.`);
        modeIndex = 1;
    }
    
    let targetEigenVectorIndex = frequenciesFEM[modeIndex - 1].index;
    let femModeShape = [0]; // Boundary condition u(0) = 0
    
    for (let i = 0; i < n_dof; i++) {
        femModeShape.push(eigenvectors[i][targetEigenVectorIndex]);
    }
    femModeShape.push(0); // Boundary condition u(L) = 0

    // Normalize FEM mode shape so maximum absolute displacement is 1
    let maxVal = Math.max(...femModeShape.map(Math.abs));
    // Ensure the first major peak is positive for easy visual comparison
    let signFactor = femModeShape.find(v => Math.abs(v) > 0.1) > 0 ? 1 : -1;
    femModeShape = femModeShape.map(val => (val / maxVal) * signFactor);

    logText += `\nExtracted and normalized mode shape for Mode ${modeIndex}.\n`;

    // 9. Generate Plotly Data
    let x_fem = [];
    for(let i = 0; i <= N; i++) x_fem.push(i * Le);

    let x_ana = [];
    let y_ana = [];
    for(let i = 0; i <= 200; i++) { // 200 points for a smooth analytical curve
        let x = (i / 200) * L;
        x_ana.push(x);
        y_ana.push(Math.sin(modeIndex * Math.PI * x / L)); 
    }

    const traceAnalytical = {
        x: x_ana,
        y: y_ana,
        mode: 'lines',
        name: `Analytical Mode ${modeIndex}`,
        line: { color: '#e74c3c', width: 3 }
    };

    const traceFEM = {
        x: x_fem,
        y: femModeShape,
        mode: 'lines+markers',
        name: `FEM Mode ${modeIndex} (N=${N})`,
        line: { color: '#2980b9', width: 2 },
        marker: { size: 8, color: '#2c3e50' }
    };

    const layout = {
        title: `Mode ${modeIndex} Shape: Displacement u(x)`,
        xaxis: { title: 'Position along bar, x (m)', zeroline: false },
        yaxis: { title: 'Normalized Displacement', zeroline: false },
        margin: { t: 50, b: 50, l: 50, r: 20 },
        hovermode: 'closest'
    };

    Plotly.newPlot('plot', [traceAnalytical, traceFEM], layout);

    // Update log
    logText += `--- CALCULATION COMPLETE ---\n`;
    document.getElementById('calcLog').innerText = logText;
}
// Run automatically on initial load
window.onload = runAnalysis;