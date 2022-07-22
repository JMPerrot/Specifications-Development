# Specifications-Development-
In this repository, we will submit programs to establish specifications for the development of leaf-clip optical sensore dicated to the monitoring of leaf biochemical content.
The PROSPECT Model is used with ANGERS database to simulate different configurations to find an optimal set of optical parameters.
## Reference : Initial situation and goal
### Initial situation
spectral range : 

-[400-900] nm for chlorophyll (CHL) and carotenoid (CAR)

-[1300-2400] nm for equivalent water thickness (EWT) and leaf mass per area (LMA)
                 
<img src="https://user-images.githubusercontent.com/101126884/180430522-05cb82e2-94ec-4c9c-839a-479ea0ef6d7a.png" width="250" height="250">

### Approximative goal
spectral range : 

[700-720] nm for chlorophyll (CHL)


[-]carotenoid (CAR)

[1300-2400] nm for equivalent water thickness (EWT) and leaf mass per area (LMA)
                 
<img src="https://user-images.githubusercontent.com/101126884/180439261-c12be332-4d4c-4ff8-bd70-d7b8dfc2f078.png" width="250" height="250">          

We can observe that we have an improvement in the normalized root mean square error between these two situations. The objective of this repository is to share codes retracing the steps allowing, from the initial situation, to tend towards the optimized situation by using specific wavelengths and not extended spectral ranges
## R&T

### First step : Sampling stepping

<img src="https://user-images.githubusercontent.com/101126884/180441664-f2309ede-4946-4674-ae48-56db45285244.png" width="250" height="250">   

### Second step : Sampling translation

<img src="https://user-images.githubusercontent.com/101126884/180441864-b0ab02b3-f158-411d-8f45-788fb4adbcf4.png" width="250" height="250">   

### Third step : features selection

<img src="https://user-images.githubusercontent.com/101126884/180441928-c20c6130-5cbc-4b90-b72a-9c06693990e9.png" width="250" height="250">   

### Fourth step : gaussian filter application

<img src="https://user-images.githubusercontent.com/101126884/180441997-f5215297-6bc1-4db8-94e2-a4af1d866f81.png" width="250" height="250">   

## R only with prior N estimation
N is the number of leaf structure in the plate leaf model

### First step : Sampling stepping

<img src="https://user-images.githubusercontent.com/101126884/180442274-12d6bdc0-411a-4f37-9704-80cf6f0ba9ae.png" width="250" height="250">   

### Second step : Sampling translation

<img src="https://user-images.githubusercontent.com/101126884/180442341-9f528c50-f31a-404e-a3dd-8e9cb45e3681.png" width="250" height="250">   

### Third step : features selection

<img src="https://user-images.githubusercontent.com/101126884/180442416-d2c7cbec-b9fe-46f8-90e0-d4a5ec4eb3b5.png" width="250" height="250">   

### Fourth step : gaussian filter application

<img src="https://user-images.githubusercontent.com/101126884/180446348-eadbbd4a-8110-4769-88e6-2405ebcd9b79.png" width="250" height="250">   

## T only with prior N estimation

### First step : Sampling stepping

<img src="https://user-images.githubusercontent.com/101126884/180442547-c51dc224-59d3-4f19-abe2-57a5846c1b98.png" width="250" height="250">   

### Second step : Sampling translation

<img src="https://user-images.githubusercontent.com/101126884/180442604-16476080-0820-49b6-9a48-7fe7296998e0.png" width="250" height="250">   

### Third step : features selection

<img src="https://user-images.githubusercontent.com/101126884/180442674-f96753d3-a570-44ba-898d-798e92858e85.png" width="250" height="250">   

### Fourth step : gaussian filter application

<img src="" width="250" height="250">   
