// Navier-Stokes 2D stable incompressible fluid solver.
// Code implements inviscid flow. The solution is Euler-like.
// There is no thermal term as well.
//
// Performance is pretty good, but code favors compactness (~500 lines including Unity boilerplating), simplicity and readability.

using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class cfd : MonoBehaviour
{
// public
    [Space][Header("grid")]
    public int resolution = 100;
    public bool boundary = false;
    [Space][Header("emission")]
    public float injectDensity = 1;
    public float injectVelocity = 1;
    public float maxInjectVelocity = 10;
    public float radius = 0;
    [Range(1,100)] public int samples = 1;
    [Space][Header("simulation")]
    [Range(0, 1)] public float friction = 0.001f;
    [Range(-5, 5)] public float gravity = 1;
    [Range(0, 1)] public float swirl = 0;
    [Space][Header("diffusion")]
    public bool diffuseDensity = true;
    public bool diffuseVelocity = true;
    [Range(0,1)] public float densityDiffusion = 1;
    [Range(0, 1)] public float velocityDiffusion = 0.1f;
    [Range(0, 5)] public int diffusionQuality = 0;
    [Space][Header("advection")]
    public bool advectDensity = true;
    public bool advectVelocity = true;
    public float densityAdvection = 1;
    public float velocityAdvection = 1;
    [Space][Header("dissipation")]
    public float densityDissipation = 0;
    public float velocityDissipation = 0;
    [Space][Header("conservation")]
    public bool conserveMass = true;
    [Range(1, 50)] public int solverQuality = 20;
    [Space][Header("visualisation")]
    public bool displayVelocity;
    public float displayVelocityMultiplier = 1;
    public enum DisplayGrid{none, density, divergence, pressure};
    public DisplayGrid displayGrid;
    public float displayGridMultiplier = 1;

// private
    int voxelCount;
    float voxelScale, voxelScale2;
    int resolution0, resolution2 = 0;
    float[] density0, density, divergence, pressure, vorticity, absoluteVorticity, vorticityLength;
    Vector2[] velocity0, velocity, vorticityGradient, vorticityConfinement;
    Vector3 v;

    ParticleSystem.Particle[] particles_scalarField, particles_vectorField;
    ParticleSystem particleSystem_scalarField, particleSystem_vectorField;
    Renderer renderer_scalarField, renderer_vectorField;

    // (Re)calculate grid resolution and voxel size, zero-out fields,
    // update the particle systems that display data on screen.
    void Initialize() 
    {
        // Run this code at start-up only, or if the user changed grid resolution.
        if (resolution0 != resolution)
        {
            Vector3 scale = transform.localScale;
            if (scale.x == 0.0f || resolution <= 0) return;

            // Voxel resolution, count and scale.
            resolution0 = resolution;
            resolution2 = resolution+2;
            voxelScale = scale.x / resolution;
            voxelScale2 = scale.x / resolution2;
            voxelCount = resolution2 * resolution2;

            // Zero-out fields.
            density0 = density = new float[voxelCount];
            velocity0 = velocity = new Vector2[voxelCount];
            divergence = new float[voxelCount];
            pressure = new float[voxelCount];
            vorticity = new float[voxelCount];
            absoluteVorticity = new float[voxelCount];
            vorticityLength = new float[voxelCount];
            vorticityGradient = new Vector2[voxelCount];
            vorticityConfinement = new Vector2[voxelCount];

            // Initialize two particle systems displaying the voxel grid on screen (its velocity, density, divergence, pressure fields).
            int ii = 0;
            float offsetX = scale.x * 0.5f - voxelScale * 0.5f;
            float offsetY = scale.y * 0.5f - voxelScale * 0.5f;
            particles_scalarField = new ParticleSystem.Particle[resolution*resolution];
            particles_vectorField = new ParticleSystem.Particle[resolution*resolution];
            for (int i = 0; i < resolution; ++i)
            for (int j = 0; j < resolution; ++j)
            {
                v.x = i * voxelScale - offsetX;
                v.y = j * voxelScale - offsetY;
                particles_scalarField[ii].position = v;
                particles_scalarField[ii].startSize = voxelScale;
                particles_vectorField[ii].position = v;
                particles_vectorField[ii].startSize = 0.01f;
                ++ii;
            }
            // Create the particles.
            var emitParams = new ParticleSystem.EmitParams();
            emitParams.position = Vector3.zero;
            particleSystem_scalarField.Emit(emitParams, resolution*resolution);
            particleSystem_scalarField.SetParticles(particles_scalarField, resolution*resolution);
            particleSystem_vectorField.Emit(emitParams, resolution*resolution);
            particleSystem_vectorField.SetParticles(particles_vectorField, resolution*resolution);
        }
    }

    // Emit density and velocity into the fluid grid.
    // The emitter location inside the grid is linked to the cursor's screen-space coordinates.
    // LMB emits density, RMB emits velocity.
    Vector2 position0 = new Vector2();
    Vector2 bend = new Vector2();
    void Emit()
    {
        if (!Input.GetMouseButton(0) && !Input.GetMouseButton(1)) return;

        RaycastHit hit;
        if (!Physics.Raycast(Camera.main.ScreenPointToRay(Input.mousePosition), out hit)) return;

        // Emit multiple times (defined by the number of "samples") within the cursor area (defined by its "radius" parameter).
        int idx;
        float s = transform.localScale.x;
        for (int i = 0; i < samples; ++i)
        {
            if (samples == 1) v = Vector3.zero;
            else v = Quaternion.AngleAxis(UnityEngine.Random.Range(0, 360), Vector3.forward) * Vector3.right * UnityEngine.Random.Range(0, radius);
            v += hit.point;
            // Skip if the emission point is outside the voxel grid.
            if (v.x < 0 || v.x > s) return;
            if (v.y < 0 || v.y > s) return;
            // Calculate the voxel index.
            v.x = (float)Math.Ceiling(v.x / voxelScale);
            v.y = (float)Math.Ceiling(v.y / voxelScale);
            idx = (int)(v.x + v.y * resolution2 + 0.5f);

            // Add density.
            if (Input.GetMouseButton(0)) density0[idx] = injectDensity;

            // Add velocity.
            if (Input.GetMouseButtonDown(1))
            {
                position0.x = hit.point.x;
                position0.y = hit.point.y;
            }
            else if (Input.GetMouseButton(1))
            {
                // Calculate direction vector from previous and current cursor positions.
                v.x = position0.x = hit.point.x - position0.x;
                v.y = position0.y = hit.point.y - position0.y;
                // Bend it towards the mouse direction, so the result feels more natural to the user.
                position0.x += bend.x;
                position0.y += bend.y;
                velocity0[idx] += position0.normalized * injectVelocity;
                // Clamp the result, so it does not overshoot.
                if (velocity0[idx].magnitude > maxInjectVelocity)
                    velocity0[idx] = velocity0[idx].normalized * maxInjectVelocity;
                // Remember the current bend vector.
                bend.x = v.x;
                bend.y = v.y;
                // Remember the current cursor position.
                position0.x = hit.point.x;
                position0.y = hit.point.y;
            }
        }
    }

    // Generate small scale vorticles to increase detail in the simulation.
    void VorticityConfinement()
    {
        if (swirl <= 0) return; else if (swirl > 1) swirl = 1;

        int ii;
        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            vorticity[ii] = (velocity0[ii+resolution2].y-velocity0[ii-resolution2].y - velocity0[ii+1].x+velocity0[ii-1].x) * 0.5f;
            if (vorticity[ii] >= 0) absoluteVorticity[ii] = vorticity[ii];
            else absoluteVorticity[ii] = -vorticity[ii];
        }
        if (boundary)
        {
            SetScalarBoundary(ref vorticity);
            SetScalarBoundary(ref absoluteVorticity);
        }

        Vector2 v;
        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            v.x = (absoluteVorticity[ii+1] - absoluteVorticity[ii-1]) * 0.5f;
            v.y = (absoluteVorticity[ii+resolution2] - absoluteVorticity[ii-resolution2]) * 0.5f;
            vorticityGradient[ii] = v;
            vorticityLength[ii] = Mathf.Sqrt(vorticityGradient[ii].x*vorticityGradient[ii].x + vorticityGradient[ii].y*vorticityGradient[ii].y);
            if (vorticityLength[ii] < 0.01f)
                vorticityConfinement[ii] = Vector2.zero;
            else
            {
                v.x = vorticityGradient[ii].x / vorticityLength[ii];
                v.x = vorticityGradient[ii].y / vorticityLength[ii];
                vorticityConfinement[ii] = v;
            }
        }
        if (boundary) SetVectorBoundary(ref vorticityConfinement);

        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            v.x = swirl * (vorticityConfinement[ii].y * vorticity[ii]);
            v.y = swirl * (-vorticityConfinement[ii].x * vorticity[ii]);
            velocity0[ii] += v;
        }
        if (boundary) SetVectorBoundary(ref velocity0);
    }

    // Diffuse the density/velocity fields.
    void Diffuse()
    {
        // Make sure that dissipation is computed (if the related parameters are non-zero) regardless of the state of the diffuse density/velocity flags.
        if(!diffuseDensity && densityDissipation > 0) diffuseDensity = true;
        if(!diffuseVelocity && velocityDissipation > 0) diffuseVelocity = true;

        if (!diffuseDensity && !diffuseVelocity) return;

        if (densityDissipation < 0) densityDissipation = 0;
        if (velocityDissipation < 0) velocityDissipation = 0;

        int ii;
        float idd = 1 - densityDiffusion;
        float idv = 1 - velocityDiffusion;
        float a = densityDiffusion * resolution * resolution * Time.deltaTime;
        float b = velocityDiffusion * resolution * resolution * Time.deltaTime;
        for(int jj = 0; jj <= diffusionQuality; ++jj)
        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            if (diffuseDensity)
            {
                density[ii] = (density0[ii] * idd + (((density0[ii] + a * (density0[ii - 1] + density0[ii + 1] + density0[ii - resolution2] + density0[ii + resolution2])) / (1 + 4 * a)) * densityDiffusion)) - densityDissipation;
                if (density[ii] < 0) density[ii] = 0;
            }
            if (diffuseVelocity)
            {
                velocity[ii] = (velocity0[ii] * idv + (((velocity0[ii] + b * (velocity0[ii - 1] + velocity0[ii + 1] + velocity0[ii - resolution2] + velocity0[ii + resolution2])) / (1 + 4 * b)) * velocityDiffusion));
                if (velocityDissipation > 0)
                {
                    float m = velocity[ii].magnitude - velocityDissipation;
                    if (m < 0) m = 0;
                    velocity[ii] = velocity[ii].normalized * m;
                }
            }
        }
        if (diffuseDensity)
        {
            density0 = density;
            if (boundary) SetScalarBoundary(ref density0);
        }
        if (diffuseVelocity)
        {
            velocity0 = velocity;
            if (boundary) SetVectorBoundary(ref velocity0);
        }
    }

    // Advecty densitiy along the velocity field using backward advection.
    void AdvectDensity()
    {
        if (!advectDensity) return;

        int ii, i0, j0, i1, j1;
        float x, y, s1, s0, t1, t0;
        float c = Time.deltaTime * resolution * densityAdvection;
        density = new float[voxelCount];
        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            x = i - velocity0[ii].x * c;
            y = j - velocity0[ii].y * c;
            if (x < 0.5f) x = 0.5f; else if (x > resolution+0.5f) x = resolution+0.5f; i0 = (int)x; i1 = i0+1;
            if (y < 0.5f) y = 0.5f; else if (y > resolution+0.5f) y = resolution+0.5f; j0 = (int)y; j1 = j0+1;
            s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
            density[ii] = s0*(t0*density0[i0+resolution2*j0] + t1*density0[i0+resolution2*j1]) + s1*(t0*density0[i1+resolution2*j0] + t1*density0[i1+resolution2*j1]);
        }
        density0 = density;
        if (boundary) SetScalarBoundary(ref density0);
    }

    // Self-advect velocity field using backward advection.
    void AdvectVelocity()
    {
        if (!advectVelocity) return;

        if (friction <= 0) friction = 0; else if (friction > 1) friction = 1;

        int ii, i0, j0, i1, j1;
        float x, y, s1, s0, t1, t0;
        float c = Time.deltaTime * resolution * velocityAdvection;
        float f = 1 - friction;
        Vector2 g = new Vector2(0, gravity * 0.001f);
        velocity = new Vector2[voxelCount];
        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            x = i - velocity0[ii].x * c;
            y = j - velocity0[ii].y * c;
            if (x < 0.5f) x = 0.5f; else if (x > resolution+0.5f) x = resolution+0.5f; i0 = (int)x; i1 = i0+1;
            if (y < 0.5f) y = 0.5f; else if (y > resolution+0.5f) y = resolution+0.5f; j0 = (int)y; j1 = j0+1;
            s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
            velocity[ii] = (s0*(t0*velocity0[i0+resolution2*j0] + t1*velocity0[i0+resolution2*j1]) + s1*(t0*velocity0[i1+resolution2*j0] + t1*velocity0[i1+resolution2*j1])) * f + density0[ii] * g;
        }
        velocity0 = velocity;
        if (boundary) SetVectorBoundary(ref velocity0);
    }

    // Mass conservation.
    void Project()
    {
        if (!conserveMass) return;

        int ii;
        float h = 1.0f / resolution;
        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            divergence[ii] = -0.5f * h * (velocity0[ii+1].x - velocity0[ii-1].x + velocity0[ii+resolution2].y - velocity0[ii-resolution2].y);
        }
        if (boundary) SetScalarBoundary(ref divergence);

        for (int k = 0; k < solverQuality; ++k)
            for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                pressure[ii] = (divergence[ii] + pressure[ii-1] + pressure[ii+1] + pressure[ii-resolution2] + pressure[ii+resolution2]) / 4;
            }
        if (boundary) SetScalarBoundary(ref pressure);

        for (int i = 1; i <= resolution; ++i)
        for (int j = 1; j <= resolution; ++j)
        {
            ii = i + resolution2 * j;
            velocity0[ii].x -= 0.5f * (pressure[ii+1] - pressure[ii-1]) / h;
            velocity0[ii].y -= 0.5f * (pressure[ii+resolution2] - pressure[ii-resolution2]) / h;
        }
        if (boundary) SetVectorBoundary(ref velocity0);
    }

    void SetScalarBoundary(ref float[] array)
    {
        int idx = resolution2 * resolution2 - resolution2;
        int idx2;
        for (int i = 1; i <= resolution; ++i)
        {
            array[i] = array[i+resolution2]; // bottom row
            array[i+idx] = array[i+idx-resolution2]; // top row
            // left column
            idx2 = i*resolution2;
            array[idx2] = array[idx2+1];
            // right column
            idx2 += resolution2;
            array[idx2-1] = array[idx2-2];
        }

        array[0] = (array[1] + array[resolution2]) * 0.5f; // bottom-left corner
        array[idx] = (array[idx+1] + array[idx-resolution2]) * 0.5f; // top-left corner
        array[resolution2-1] = (array[resolution2-2] + array[resolution2*2-1]) * 0.5f; // bottom-right corner
        idx = resolution2*resolution2;
        array[idx-1] = (array[idx-2] + array[idx-resolution2-1])/2; //top-right corner
    }

    void SetVectorBoundary(ref Vector2[] array)
    {
        Vector2 v;
        int idx = resolution2 * resolution2 - resolution2;
        int idx2;
        for(int i = 1; i <= resolution; ++i)
        {
            // bottom row
            v.x = array[i+resolution2].x;
            v.y = -array[i+resolution2].y; 
            array[i] = v;
            // top row
            v.x = array[i+idx-resolution2].x;
            v.y = -array[i+idx-resolution2].y; 
            array[i+idx] = v;
            // left column
            idx2 = i*resolution2;
            v.x = -array[idx2+1].x;
            v.y = array[idx2+1].y; 
            array[idx2] = v;
            // right column
            idx2 += resolution2;
            v.x = -array[idx2-2].x;
            v.y = array[idx2-2].y;
            array[idx2-1] = v;
        }

        array[0] = (array[1] + array[resolution2]) * 0.5f; // bottom-left corner
        array[idx] = (array[idx+1] + array[idx-resolution2]) * 0.5f; // top-left corner
        array[resolution2-1] = (array[resolution2-2] + array[resolution2*2-1]) * 0.5f; // bottom-right corner
        idx = resolution2*resolution2;
        array[idx-1] = (array[idx-2] + array[idx-resolution2-1])/2; //top-right corner
    }

    // Draw voxel grid on screen.
    void Display()
    {
        // Scalar fields.
        if (displayGrid != DisplayGrid.none)
        {
            int ii;
            float value;
            renderer_scalarField.enabled = true;
            if (displayGrid == DisplayGrid.density)
            {
                ii = 0;
                for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    value = density[i+resolution2*j] * displayGridMultiplier;
                    particles_scalarField[ii].startColor = new Color(value, value, value, 1);
                    ++ii;
                }
            }
            else if (displayGrid == DisplayGrid.divergence)
            {
                ii = 0;
                for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    value = divergence[i+resolution2*j] * displayGridMultiplier * 10000;
                    particles_scalarField[ii].startColor = new Color(value, value, value, 1);
                    ++ii;
                }
            }
            else if (displayGrid == DisplayGrid.pressure)
            {
                ii = 0;
                for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    value = pressure[i+resolution2*j] * displayGridMultiplier * 10000;
                    particles_scalarField[ii].startColor = new Color(value, value, value, 1);
                    ++ii;
                }
            }
            particleSystem_scalarField.SetParticles(particles_scalarField, voxelCount);
        }
        else renderer_scalarField.enabled = false;

        // Vector field (in this case that's the velocity only).
        if (displayVelocity)
        {
            int ii = 0;
            renderer_vectorField.enabled = true;
            for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                particles_vectorField[ii].velocity = -velocity0[i+resolution2*j] * displayVelocityMultiplier;
                particles_vectorField[ii].startColor = new Color(1,1,1,1);
                ++ii;
            }
            particleSystem_vectorField.SetParticles(particles_vectorField, voxelCount);
        }
        else renderer_vectorField.enabled = false;
    }

    void Start()
    {
        transform.position = new Vector3(Mathf.Abs(transform.localScale.x)*0.5f, Mathf.Abs(transform.localScale.x)*0.5f, 0);
        transform.localScale = new Vector3(Mathf.Abs(transform.localScale.x), Mathf.Abs(transform.localScale.x), 1);
        particleSystem_scalarField = transform.GetChild(0).GetComponent<ParticleSystem>();
        renderer_scalarField = particleSystem_scalarField.GetComponent<Renderer>();
        particleSystem_vectorField = transform.GetChild(1).GetComponent<ParticleSystem>();
        renderer_vectorField = particleSystem_vectorField.GetComponent<Renderer>();
    }

    void Update()
    {
        Initialize();
        Emit();
        Diffuse();
        Project();
        AdvectVelocity();
        Project();
        AdvectDensity();
        VorticityConfinement();
        Display();
    }
}