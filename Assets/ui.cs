// Custom-tailored UI to the cfd solver's public properties.
// Code is based on the hierarchy order of the "ui" prefab.

using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ui : MonoBehaviour
{
    public cfd solver;
    [Header("ranges")]
    public Vector2 resolution = new Vector2(10, 400);
    public Vector2 radius = new Vector2(0, 10);
    public Vector2 samples = new Vector2(0, 100);
    public Vector2 friction = new Vector2(0, 1);
    public Vector2 gravity = new Vector2(-5, 5);
    public Vector2 swirl = new Vector2(0, 1);
    public Vector2 densityDiffusion = new Vector2(0, 1);
    public Vector2 velocityDiffusion = new Vector2(0, 1);
    public Vector2 diffusionQuality = new Vector2(0, 10);
    public Vector2 solverQuality = new Vector2(1, 50);

    Text framerate;
    float deltaTime = 0.0f;
    float refreshInterval = 0.25f;
    float time = 0;

    public void Start()
    {
        framerate = transform.GetChild(0).GetChild(transform.GetChild(0).childCount - 1).GetChild(0).GetComponent<Text>();

        Slider slider;
        
        transform.GetChild(0).GetChild(0).GetChild(1).GetComponent<InputField>().text = solver.resolution.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(2).GetComponent<Slider>();
        slider.value = solver.resolution;
        slider.minValue = resolution.x;
        slider.maxValue = resolution.y;

        transform.GetChild(0).GetChild(1).GetComponent<Toggle>().isOn = solver.boundary;

        transform.GetChild(0).GetChild(2).GetChild(1).GetComponent<InputField>().text = solver.injectDensity.ToString();

        transform.GetChild(0).GetChild(3).GetChild(1).GetComponent<InputField>().text = solver.injectVelocity.ToString();

        transform.GetChild(0).GetChild(4).GetChild(1).GetComponent<InputField>().text = solver.maxInjectVelocity.ToString();

        transform.GetChild(0).GetChild(5).GetChild(1).GetComponent<InputField>().text = solver.radius.ToString();
        slider = transform.GetChild(0).GetChild(5).GetChild(2).GetComponent<Slider>();
        slider.value = solver.radius;
        slider.minValue = radius.x;
        slider.maxValue = radius.y;

        transform.GetChild(0).GetChild(6).GetChild(1).GetComponent<InputField>().text = solver.samples.ToString();
        slider = transform.GetChild(0).GetChild(6).GetChild(2).GetComponent<Slider>();
        slider.value = solver.samples;
        slider.minValue = samples.x;
        slider.maxValue = samples.y;

        transform.GetChild(0).GetChild(7).GetChild(1).GetComponent<InputField>().text = solver.friction.ToString();
        slider = transform.GetChild(0).GetChild(7).GetChild(2).GetComponent<Slider>();
        slider.value = solver.friction;
        slider.minValue = friction.x;
        slider.maxValue = friction.y;

        transform.GetChild(0).GetChild(8).GetChild(1).GetComponent<InputField>().text = solver.gravity.ToString();
        slider = transform.GetChild(0).GetChild(8).GetChild(2).GetComponent<Slider>();
        slider.value = solver.gravity;
        slider.minValue = gravity.x;
        slider.maxValue = gravity.y;

        transform.GetChild(0).GetChild(9).GetChild(1).GetComponent<InputField>().text = solver.swirl.ToString();
        slider = transform.GetChild(0).GetChild(9).GetChild(2).GetComponent<Slider>();
        slider.value = solver.swirl;
        slider.minValue = swirl.x;
        slider.maxValue = swirl.y;

        transform.GetChild(0).GetChild(10).GetComponent<Toggle>().isOn = solver.diffuseDensity;
        transform.GetChild(0).GetChild(11).GetComponent<Toggle>().isOn = solver.diffuseVelocity;

        transform.GetChild(0).GetChild(12).GetChild(1).GetComponent<InputField>().text = solver.densityDiffusion.ToString();
        slider = transform.GetChild(0).GetChild(12).GetChild(2).GetComponent<Slider>();
        slider.value = solver.densityDiffusion;
        slider.minValue = densityDiffusion.x;
        slider.maxValue = densityDiffusion.y;

        transform.GetChild(0).GetChild(13).GetChild(1).GetComponent<InputField>().text = solver.velocityDiffusion.ToString();
        slider = transform.GetChild(0).GetChild(13).GetChild(2).GetComponent<Slider>();
        slider.value = solver.velocityDiffusion;
        slider.minValue = velocityDiffusion.x;
        slider.maxValue = velocityDiffusion.y;

        transform.GetChild(0).GetChild(14).GetChild(1).GetComponent<InputField>().text = solver.diffusionQuality.ToString();
        slider = transform.GetChild(0).GetChild(14).GetChild(2).GetComponent<Slider>();
        slider.minValue = diffusionQuality.x;
        slider.maxValue = diffusionQuality.y;
        slider.value = solver.diffusionQuality;

        transform.GetChild(0).GetChild(15).GetComponent<Toggle>().isOn = solver.advectDensity;
        transform.GetChild(0).GetChild(16).GetComponent<Toggle>().isOn = solver.advectVelocity;

        transform.GetChild(0).GetChild(17).GetChild(1).GetComponent<InputField>().text = solver.densityAdvection.ToString();
        transform.GetChild(0).GetChild(18).GetChild(1).GetComponent<InputField>().text = solver.velocityAdvection.ToString();
        transform.GetChild(0).GetChild(19).GetChild(1).GetComponent<InputField>().text = solver.densityDissipation.ToString();
        transform.GetChild(0).GetChild(20).GetChild(1).GetComponent<InputField>().text = solver.velocityDissipation.ToString();

        transform.GetChild(0).GetChild(21).GetComponent<Toggle>().isOn = solver.conserveMass;

        transform.GetChild(0).GetChild(22).GetChild(1).GetComponent<InputField>().text = solver.solverQuality.ToString();
        slider = transform.GetChild(0).GetChild(22).GetChild(2).GetComponent<Slider>();
        slider.value = solver.solverQuality;
        slider.minValue = solverQuality.x;
        slider.maxValue = solverQuality.y;
 
        transform.GetChild(0).GetChild(23).GetComponent<Toggle>().isOn = solver.displayVelocity;

        transform.GetChild(0).GetChild(24).GetChild(1).GetComponent<InputField>().text = solver.displayVelocityMultiplier.ToString();

        transform.GetChild(0).GetChild(25).GetComponent<Dropdown>().value = (int)(solver.displayGrid);

        transform.GetChild(0).GetChild(26).GetChild(1).GetComponent<InputField>().text = solver.displayGridMultiplier.ToString();
    }

    void Update()
    {
        if (Input.GetKeyDown("p")) solver.enabled = !solver.enabled;

        // Display framerate.
        deltaTime += (Time.unscaledDeltaTime - deltaTime) * 0.1f;
        time += Time.unscaledDeltaTime;
        if (time < refreshInterval) return;
        time = 0;
        float msec = deltaTime * 1000.0f;
        float fps = 1.0f / deltaTime;
        if (fps > 28) framerate.color = Color.green;
        else if (fps >= 25 && fps < 28) framerate.color = Color.yellow;
        else framerate.color = Color.red;
        framerate.text = string.Format("framerate {1:0.} fps ({0:0.0} ms)", msec, fps);
    }

    void UpdateSlider(InputField field)
    {
        field.transform.parent.GetChild(2).GetComponent<Slider>().value = Convert.ToSingle(field.text);
    }

    void UpdateField(Slider slider)
    {
        slider.transform.parent.GetChild(1).GetComponent<InputField>().text = ((float)(slider.value)).ToString();
    }

    public void setResolutionFromSlider(Slider slider)
    {
        solver.resolution = (int)slider.value;
        UpdateField(slider);
    }
    public void setResolutionFromField(InputField field)
    {
        solver.resolution = (int)Mathf.Clamp(Convert.ToInt32(field.text), resolution.x, resolution.y);
        field.text = solver.resolution.ToString();
        UpdateSlider(field);
    }

    public void setBoundary(Toggle toggle)
    {
        solver.boundary = toggle.isOn;
    }

    public void setDensity(InputField field)
    {
        solver.injectDensity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.injectDensity.ToString();
    }

    public void setVelocity(InputField field)
    {
        solver.injectVelocity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.injectVelocity.ToString();
    }

    public void setMaxVelocity(InputField field)
    {
        solver.maxInjectVelocity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.maxInjectVelocity.ToString();
    }

    public void setRadiusFromSlider(Slider slider)
    {
        solver.radius = (float)slider.value;
        UpdateField(slider);
    }

    public void setRadiusFromField(InputField field)
    {
        solver.radius = Mathf.Clamp(Convert.ToSingle(field.text), radius.x, radius.y);
        field.text = solver.radius.ToString();
        UpdateSlider(field);
    }

    public void setSamplesFromSlider(Slider slider)
    {
        solver.samples = (int)slider.value;
        UpdateField(slider);
    }
    public void setSamplesFromField(InputField field)
    {
        solver.samples = (int)Mathf.Clamp(Convert.ToInt32(field.text), samples.x, samples.y);
        field.text = solver.samples.ToString();
        UpdateSlider(field);
    }

    public void setFrictionFromSlider(Slider slider)
    {
        solver.friction = slider.value;
        UpdateField(slider);
    }
    public void setFrictionFromField(InputField field)
    {
        solver.friction = (float)Mathf.Clamp(Convert.ToSingle(field.text), friction.x, friction.y);
        field.text = solver.friction.ToString();
        UpdateSlider(field);
    }

    public void setGravityFromSlider(Slider slider)
    {
        solver.gravity = slider.value;
        UpdateField(slider);
    }

    public void setGravityFromField(InputField field)
    {
        solver.gravity = (float)Mathf.Clamp(Convert.ToSingle(field.text), gravity.x, gravity.y);
        field.text = solver.gravity.ToString();
        UpdateSlider(field);
    }

    public void setSwirlFromSlider(Slider slider)
    {
        solver.swirl = slider.value;
        UpdateField(slider);
    }

    public void setSwirlFromField(InputField field)
    {
        solver.swirl = (float)Mathf.Clamp(Convert.ToSingle(field.text), swirl.x, swirl.y);
        field.text = solver.swirl.ToString();
        UpdateSlider(field);
    }

    public void diffuseDensity(Toggle toggle)
    {
        solver.diffuseDensity = toggle.isOn;
    }

    public void diffuseVelocity(Toggle toggle)
    {
        solver.diffuseVelocity = toggle.isOn;
    }

    public void setDensityDiffusionFromSlider(Slider slider)
    {
        solver.densityDiffusion = slider.value;
        UpdateField(slider);
    }

    public void setDensityDiffusionFromField(InputField field)
    {
        solver.densityDiffusion = (float)Mathf.Clamp(Convert.ToSingle(field.text), densityDiffusion.x, densityDiffusion.y);
        field.text = solver.densityDiffusion.ToString();
        UpdateSlider(field);
    }

    public void setVelocityDiffusionFromSlider(Slider slider)
    {
        solver.velocityDiffusion = slider.value;
        UpdateField(slider);
    }

    public void setVelocityDiffusionFromField(InputField field)
    {
        solver.velocityDiffusion = (float)Mathf.Clamp(Convert.ToSingle(field.text), velocityDiffusion.x, velocityDiffusion.y);
        field.text = solver.velocityDiffusion.ToString();
        UpdateSlider(field);
    }

    public void setDiffusionQualityFromSlider(Slider slider)
    {
        solver.diffusionQuality = (int)slider.value;
        UpdateField(slider);
    }

    public void setDiffusionQualityFromField(InputField field)
    {
        solver.diffusionQuality = (int)Mathf.Clamp(Convert.ToInt32(field.text), diffusionQuality.x, diffusionQuality.y);
        field.text = solver.diffusionQuality.ToString();
        UpdateSlider(field);
    }

    public void advectDensity(Toggle toggle)
    {
        solver.advectDensity = toggle.isOn;
    }

    public void advectVelocity(Toggle toggle)
    {
        solver.advectVelocity = toggle.isOn;
    }

    public void setDensityAdvection(InputField field)
    {
        solver.densityAdvection = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.densityAdvection.ToString();
    }

    public void setVelocityAdvection(InputField field)
    {
        solver.velocityAdvection = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.velocityAdvection.ToString();
    }

    public void setDensityDissipation(InputField field)
    {
        solver.densityDissipation = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.densityDissipation.ToString();
    }

    public void setVelocityDissipation(InputField field)
    {
        solver.velocityDissipation = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.velocityDissipation.ToString();
    }

    public void conserveMass(Toggle toggle)
    {
        solver.conserveMass = toggle.isOn;
    }

    public void setSolverQualityFromSlider(Slider slider)
    {
        solver.solverQuality = (int)slider.value;
        UpdateField(slider);
    }

    public void setSolverQualityFromField(InputField field)
    {
        solver.solverQuality = (int)Mathf.Clamp(Convert.ToInt32(field.text), solverQuality.x, solverQuality.y);
        field.text = solver.solverQuality.ToString();
        UpdateSlider(field);
    }

    public void displayVelocity(Toggle toggle)
    {
        solver.displayVelocity = toggle.isOn;
    }

    public void setVelocityMultiplier(InputField field)
    {
        solver.displayVelocityMultiplier = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.displayVelocityMultiplier.ToString();
    }

    public void setGridDisplay(Dropdown dropdown)
    {
        if (dropdown.value == 0) solver.displayGrid = cfd.DisplayGrid.none;
        else if (dropdown.value == 1) solver.displayGrid = cfd.DisplayGrid.density;
        else if (dropdown.value == 2) solver.displayGrid = cfd.DisplayGrid.divergence;
        else if (dropdown.value == 3) solver.displayGrid = cfd.DisplayGrid.pressure;
    }

    public void setDisplayGridMultiplier(InputField field)
    {
        solver.displayGridMultiplier = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = solver.displayGridMultiplier.ToString();
    }
}
