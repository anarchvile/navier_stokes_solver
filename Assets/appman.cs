using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;

public class appman : MonoBehaviour
{
    int scene;

    void Awake()
    {
        // Limit framerate to 30 fps.
        QualitySettings.vSyncCount = 0;
        Application.targetFrameRate = 30;
    }

    void Start ()
    {
        // Get available scenes.
        int c = SceneManager.sceneCountInBuildSettings;
        List<string> scenes = new List<string>(c);
        for (int i = 1; i < c; ++i) scenes.Add(SceneNameFromPath(SceneUtility.GetScenePathByBuildIndex(i)));
        // Update the drop-down menut with them.
        FindObjectOfType<Dropdown>().AddOptions(scenes);

        // Load the first scene by default.
        scene = 1;
        SceneManager.LoadScene(1, LoadSceneMode.Additive);
    }

    string SceneNameFromPath(string path)
    {
        var start = path.LastIndexOf("/", StringComparison.Ordinal) + 1;
        var end = path.LastIndexOf(".", StringComparison.Ordinal);
        return path.Substring(start, end-start);
    }

    // Invoked by the "presets" drop-down menu, when user selects another scene from the list.
    public void LoadScene (Dropdown dropdown)
    {
        SceneManager.UnloadSceneAsync(scene);
        scene = dropdown.value + 1;
        SceneManager.LoadScene(scene, LoadSceneMode.Additive);
    }

    void Update()
    {
        if (Input.GetKeyDown(KeyCode.F5))
        {
            SceneManager.UnloadSceneAsync(scene);
            SceneManager.LoadScene(scene, LoadSceneMode.Additive);
        }
    }
}
