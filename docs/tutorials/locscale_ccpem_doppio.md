# Running LocScale via CCPEM Doppio

```LocScale 2``` can be run via the [CCPEM Doppio](https://www.ccpem.ac.uk/docs/software/) interface. CCPEM Doppio support all four ```LocScale 2 modes``

### Step-by-step tutorial

In this example we will use the _Rattus norvegicus_ TRPV1 channel [EMDB 5778](https://www.ebi.ac.uk/emdb/EMD-5778) from the [EMDB map
and model challenge](https://zenodo.org/records/1185426). You can [download the tutorial files here](https://surfdrive.surf.nl/files/index.php/s/PkCe60os0Vc0HYh/download).

<div class="grid cards" markdown>

-   :material-numeric-1-box:{ .lg .top } __Open Doppio GUI__  

    ---
    ![Doppio_tutorial_01](img/LocScale_tutorial_01.png)
    After starting up the [CCPEM Doppio](https://www.ccpem.ac.uk/docs/software/) interface, start a new project.

-   :material-numeric-2-box:{ .lg .top } __Locate locScale job node__  

    ---
    ![Doppio_tutorial_02](img/LocScale_tutorial_02.png)
    ```LocScale``` is located under the _Map postprocessing_ tab in the program menu.

-   :material-numeric-3-box:{ .lg .top } __LocScale GUI__  
  
    ---
    ![Doppio_tutorial_03](img/LocScale_tutorial_03.png)
    The LocScale GUI provides dynamic access to set all relevant parameters depending on the ```LocScale``` mode used.


-   :material-numeric-4-box:{ .lg .top } __Advanced options ON__  

    ---
    ![Doppio_tutorial_04](img/LocScale_tutorial_04.png)
    You can toggle ON/OFF the Advanced Options menu. In most cases the default options will be just fine.

-   :material-numeric-5-box:{ .lg .top } __Advanced options OFF__  

    ---
    ![Doppio_tutorial_05](img/LocScale_tutorial_05.png)
    We recommend to only change default parameters if really necessary. Here we will keep them turned off.

-   :material-numeric-6-box:{ .lg .top } __Select LocScale mode__  

    ---
    ![Doppio_tutorial_06](img/LocScale_tutorial_06.png)

-   :material-numeric-7-box:{ .lg .top } __Upload all files__  
    This is hopw it starts

    ---
    ![Doppio_tutorial_07](img/LocScale_tutorial_07.png)

-   :material-numeric-8-box:{ .lg .top } __Set number of CPUs__  
    This is how it continues 

    ---
    ![Doppio_tutorial_08](img/LocScale_tutorial_08.png)

-   :material-numeric-9-box:{ .lg .top } __Run LocScalle__  
    This is hopw it starts

    ---
    ![Doppio_tutorial_09](img/LocScale_tutorial_09.png)

-   :material-numeric-1-box:{ .lg .top }:material-numeric-0-box:{ .lg .top } __Inspect Results__  
    This is how it continues 

    ---
    ![Doppio_tutorial_10](img/LocScale_tutorial_10.png)

-   :material-numeric-1-box:{ .lg .top }:material-numeric-1-box:{ .lg .top } __XXX__  
    View the slices of the input half-maps and the LocScale output.

    ---
    ![Doppio_tutorial_11](img/LocScale_tutorial_11.png)

-   :material-numeric-1-box:{ .lg .top }:material-numeric-2-box:{ .lg .top } __YYY__  
    Visualise the 3D structure of the output. Change the isosurface threshold for better viewing. Typical range for maps produced by         LocScale for best visualisation is between 0.05 and 0.15

    ---
    ![Doppio_tutorial_12](img/LocScale_tutorial_12.png)

</div>

#### Relevant advanced options in CCPEM Doppio
| Option                  | Notes                                                            | Affected Method |
| ------------------------| ------------------------------------ |---------------------------|
| `LocScale Window size`  | Choose an even number, preferable range is between 15 and 30     | ALL (:material-fruit-pineapple:, :material-fruit-watermelon:, :material-fruit-cherry:, :material-fruit-pear: |
| `Resolution`            | :material-check-all: Update resource | :material-fruit-watermelon: |
| `DELETE`                | :material-close:     Delete resource | :material-fruit-watermelon: |
