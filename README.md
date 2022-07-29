# Specifications-Development-
In this repository, we will submit programs to establish specifications for the development of leaf-clip optical sensore dicated to the monitoring of leaf biochemical content.
The PROSPECT Model is used with ANGERS database to simulate different configurations to find an optimal set of optical parameters.
## Reference : Initial situation and goal
### Initial situation
spectral range : 

-[400-900] nm for chlorophyll (CHL) and carotenoid (CAR)

-[1300-2400] nm for equivalent water thickness (EWT) and leaf mass per area (LMA)
                 
<img src="https://user-images.githubusercontent.com/101126884/181711347-724ac998-6576-4862-adb7-9f3bdea6c1b5.png" width="250" height="250">

In this case we have the situations from which we will start. We will seek to improve the estimates before degrading the spectral resolution
### Approximative goal
The approximate objective is defined by the optimized estimation of the PROSPECT inversion. This inversion uses the following spectral ranges for the estimation of the different constituents : 


[700-720] nm for chlorophyll (CHL)


[520-560] nm carotenoid (CAR)


[1700-2400] nm for equivalent water thickness (EWT) and leaf mass per area (LMA)
                 
<img src="https://user-images.githubusercontent.com/101126884/180439261-c12be332-4d4c-4ff8-bd70-d7b8dfc2f078.png" width="250" height="250">          

We can observe that we have an improvement in the normalized root mean square error between these two situations. The objective of this repository is to share codes retracing the steps allowing, from the initial situation, to tend towards the optimized situation by using specific wavelengths and not extended spectral ranges
## R&T

### First step : Sampling stepping

the first step is to downsample the spectral ranges by gradually increasing the step. The NRMSE will be calculated to quantify the quality of the estimates.

<img src="https://user-images.githubusercontent.com/101126884/180441664-f2309ede-4946-4674-ae48-56db45285244.png" width="250" height="250">   

### Second step : Sampling translation

The second step will be to drag the sampling comb over the spectral range using the step minimizing the NRMSE in the first step.
The goal is to see if there is more relevant information between two wavelengths than that obtained in the first step.

<img src="https://user-images.githubusercontent.com/101126884/180441864-b0ab02b3-f158-411d-8f45-788fb4adbcf4.png" width="250" height="250">   

### Third step : features selection

The third step consists in carrying out a selection among the wavelengths constituting the sampling comb in the second step. We will therefore start from all the wavelengths (N wavelengths), delete one, calculate the NRMSE, put it back and delete another. We will keep the batch of (N-1) wavelengths and start again until there is only one wavelength left

<img src="https://user-images.githubusercontent.com/101126884/180441928-c20c6130-5cbc-4b90-b72a-9c06693990e9.png" width="250" height="250">   

### Fourth step : gaussian filter application

In the fourth step, we will apply a Gaussian to the batches of wavelengths identified in step 3. This will make it possible to simulate a source of more or less good resolution. In addition, this will allow us to conclude as to the resolution necessary to have an estimate of acceptable leaf constitution.

<img src="https://user-images.githubusercontent.com/101126884/180789867-77e61247-630a-4474-9a6c-c8dafbb38d9e.png" width="250" height="250">   

## R only with prior N estimation
N is the number of leaf structure in the plate leaf model

### First step : Sampling stepping

<img src="https://user-images.githubusercontent.com/101126884/180442274-12d6bdc0-411a-4f37-9704-80cf6f0ba9ae.png" width="250" height="250">   

### Second step : Sampling translation

<img src="https://user-images.githubusercontent.com/101126884/180442341-9f528c50-f31a-404e-a3dd-8e9cb45e3681.png" width="250" height="250">   

### Third step : features selection

<img src="https://user-images.githubusercontent.com/101126884/180442416-d2c7cbec-b9fe-46f8-90e0-d4a5ec4eb3b5.png" width="250" height="250">   

### Fourth step : gaussian filter application

<img src="https://user-images.githubusercontent.com/101126884/181703139-9193bece-da56-4080-b9b6-7ca9e4d291de.png" width="250" height="250">   

## T only with prior N estimation

### First step : Sampling stepping

<img src="https://user-images.githubusercontent.com/101126884/180442547-c51dc224-59d3-4f19-abe2-57a5846c1b98.png" width="250" height="250">   

### Second step : Sampling translation

<img src="https://user-images.githubusercontent.com/101126884/180442604-16476080-0820-49b6-9a48-7fe7296998e0.png" width="250" height="250">   

### Third step : features selection

<img src="https://user-images.githubusercontent.com/101126884/180442674-f96753d3-a570-44ba-898d-798e92858e85.png" width="250" height="250">   

### Fourth step : gaussian filter application

<img src="https://user-images.githubusercontent.com/101126884/181219598-13aff664-4730-4028-b530-427e3a50bced.png" width="250" height="250"> 

