# Function to calculate BMI_z based on gender and age
calculate_zbmi <- function(gender, age, bmi) {
  
  if (is.na(gender) | is.na(age) | is.na(bmi) ){
    output <- NA
  } else {
    # For gender 1 (Male)
    if (gender == 1) {
      
      if (age >= 2 & age < 2.5) {
        output <- ((bmi / 16.482)^(-0.624) - 1) / (-0.624 * 0.07950)
      } else if (age >= 2.5 & age < 3) {
        output <- ((bmi / 16.237)^(-0.758) - 1) / (-0.758 * 0.07911)
      } else if (age >= 3 & age < 3.5) {
        output <- ((bmi / 16.019)^(-0.888) - 1) / (-0.888 * 0.07892)
      } else if (age >= 3.5 & age < 4) {
        output <- ((bmi / 15.831)^(-1.012) - 1) / (-1.012 * 0.07905)
      } else if (age >= 4 & age < 4.5) {
        output <- ((bmi / 15.676)^(-1.130) - 1) / (-1.130 * 0.07968)
      } else if (age >= 4.5 & age < 5) {
        output <- ((bmi / 15.550)^(-1.240) - 1) / (-1.240 * 0.08107)
      } else if (age >= 5 & age < 5.5) {
        output <- ((bmi / 15.452)^(-1.342) - 1) / (-1.342 * 0.08353)
      } else if (age >= 5.5 & age < 6) {
        output <- ((bmi / 15.378)^(-1.436) - 1) / (-1.436 * 0.08729)
      } else if (age >= 6 & age < 6.5) {
        output <- ((bmi / 15.336)^(-1.517) - 1) / (-1.517 * 0.09195)
      } else if (age >= 6.5 & age < 7) {
        output <- ((bmi / 15.338)^(-1.581) - 1) / (-1.581 * 0.09685)
      } else if (age >= 7 & age < 7.5) {
        output <- ((bmi / 15.392)^(-1.621) - 1) / (-1.621 * 0.10136)
      } else if (age >= 7.5 & age < 8) {
        output <- ((bmi / 15.498)^(-1.639) - 1) / (-1.639 * 0.10532)
      } else if (age >= 8 & age < 8.5) {
        output <- ((bmi / 15.643)^(-1.644) - 1) / (-1.644 * 0.10899)
      } else if (age >= 8.5 & age < 9) {
        output <- ((bmi / 15.812)^(-1.643) - 1) / (-1.643 * 0.11270)
      } else if (age >= 9 & age < 9.5) {
        output <- ((bmi / 15.996)^(-1.641) - 1) / (-1.641 * 0.11654)
      } else if (age >= 9.5 & age < 10) {
        output <- ((bmi / 16.192)^(-1.636) - 1) / (-1.636 * 0.12030)
      } else if (age >= 10 & age < 10.5) {
        output <- ((bmi / 16.400)^(-1.625) - 1) / (-1.625 * 0.12372)
      } else if (age >= 10.5 & age < 11) {
        output <- ((bmi / 16.617)^(-1.608) - 1) / (-1.608 * 0.12660)
      } else if (age >= 11 & age < 11.5) {
        output <- ((bmi / 16.846)^(-1.586) - 1) / (-1.586 * 0.12888)
      } else if (age >= 11.5 & age < 12) {
        output <- ((bmi / 17.087)^(-1.565) - 1) / (-1.565 * 0.13057)
      } else if (age >= 12 & age < 12.5) {
        output <- ((bmi / 17.342)^(-1.547) - 1) / (-1.547 * 0.13169)
      } else if (age >= 12.5 & age < 13) {
        output <- ((bmi / 17.612)^(-1.534) - 1) / (-1.534 * 0.13227)
      } else if (age >= 13 & age < 13.5) {
        output <- ((bmi / 17.898)^(-1.524) - 1) / (-1.524 * 0.13239)
      } else if (age >= 13.5 & age < 14) {
        output <- ((bmi / 18.198)^(-1.516) - 1) / (-1.516 * 0.13211)
      } else if (age >= 14 & age < 14.5) {
        output <- ((bmi / 18.511)^(-1.509) - 1) / (-1.509 * 0.13151)
      } else if (age >= 14.5 & age < 15) {
        output <- ((bmi / 18.827)^(-1.504) - 1) / (-1.504 * 0.13064)
      } else if (age >= 15 & age < 15.5) {
        output <- ((bmi / 19.138)^(-1.498) - 1) / (-1.498 * 0.12958)
      } else if (age >= 15.5 & age < 16) {
        output <- ((bmi / 19.437)^(-1.494) - 1) / (-1.494 * 0.12839)
      } else if (age >= 16 & age < 16.5) {
        output <- ((bmi / 19.725)^(-1.490) - 1) / (-1.490 * 0.12720)
      } else if (age >= 16.5 & age < 17) {
        output <- ((bmi / 20.001)^(-1.487) - 1) / (-1.487 * 0.12610)
      } else if (age >= 17 & age < 17.5) {
        output <- ((bmi / 20.265)^(-1.486) - 1) / (-1.486 * 0.12519)
      } else if (age >= 17.5 & age < 18) {
        output <- ((bmi / 20.518)^(-1.486) - 1) / (-1.486 * 0.12448)
      } else if (age >= 18) {
        output <- ((bmi / 20.759)^(-1.487) - 1) / (-1.487 * 0.12395)
      }# End if
      
    } else if (gender==2){
      
      if (age >=2 & age <2.5) {
        output <- ((bmi/16.206)^(-0.816)-1)/(-0.816*0.08447) 
      } else if (age >=2.5 & age <3){
        output <- ((bmi/15.983)^(-0.928)-1)/(-0.928*0.08417)
      } else if (age >=3 & age <3.5){
        output <- ((bmi/15.793)^(-1.029)-1)/(-1.029*0.08424) 
      } else if(age >=3.5 & age <4){
        output <- ((bmi/15.628)^(-1.121)-1)/(-1.121*0.08476) 
      } else if(age >=4 & age <4.5){
        output <- ((bmi/15.481)^(-1.207)-1)/(-1.207*0.08580) 
      } else if(age >=4.5 & age <5){
        output <- ((bmi/15.356)^(-1.286)-1)/(-1.286*0.08755) 
      } else if (age >=5 & age <5.5){
        output <- ((bmi/15.255)^(-1.356)-1)/(-1.356*0.09019) 
      } else if(age >=5.5 & age <6){
        output <- ((bmi/15.183)^(-1.414)-1)/(-1.414*0.09383) 
      } else if(age >=6 & age <6.5){
        output <- ((bmi/15.148)^(-1.461)-1)/(-1.461*0.09820) 
      } else if(age >=6.5 & age <7){
        output <- ((bmi/15.162)^(-1.498)-1)/(-1.498*0.10286) 
      } else if(age >=7 & age <7.5){
        output <- ((bmi/15.236)^(-1.526)-1)/(-1.526*0.10735)
      } else if(age >=7.5 & age <8){
        output <- ((bmi/15.366)^(-1.546)-1)/(-1.546*0.11154)
      } else if(age >=8 & age <8.5){
        output <- ((bmi/15.534)^(-1.559)-1)/(-1.559*0.11560)
      } else if(age >=8.5 & age <9){
        output <- ((bmi/15.723)^(-1.563)-1)/(-1.563*0.11969)
      } else if(age >=9 & age <9.5){
        output <- ((bmi/15.923)^(-1.560)-1)/(-1.560*0.12384)
      } else if(age >=9.5 & age <10){
        output <- ((bmi/16.142)^(-1.549)-1)/(-1.549*0.12786)
      } else if(age >=10 & age <10.5){
        output <- ((bmi/16.387)^(-1.531)-1)/(-1.531*0.13150)
      } else if(age >=10.5 & age <11){
        output <- ((bmi/16.665)^(-1.507)-1)/(-1.507*0.13457)
      } else if(age >=11 & age <11.5){
        output <- ((bmi/16.974)^(-1.481)-1)/(-1.481*0.13702)
      } else if(age >=11.5 & age <12){
        output <- ((bmi/17.309)^(-1.456)-1)/(-1.456*0.13884)
      } else if(age >=12 & age <12.5){
        output <- ((bmi/17.665)^(-1.435)-1)/(-1.435*0.14002)
      } else if(age >=12.5 & age <13){
        output <- ((bmi/18.033)^(-1.421)-1)/(-1.421*0.14057)
      } else if(age >=13 & age <13.5){
        output <- ((bmi/18.400)^(-1.412)-1)/(-1.412*0.14055)
      } else if(age >=13.5 & age <14){
        output <- ((bmi/18.755)^(-1.407)-1)/(-1.407*0.13999)
      } else if(age >=14 & age <14.5){
        output <- ((bmi/19.090)^(-1.408)-1)/(-1.408*0.13897)
      } else if(age >=14.5 & age <15){
        output <- ((bmi/19.401)^(-1.412)-1)/(-1.412*0.13763)
      } else if(age >=15 & age <15.5){
        output <- ((bmi/19.683)^(-1.419)-1)/(-1.419*0.13609)
      } else if(age >=15.5 & age <16){
        output <- ((bmi/19.934)^(-1.427)-1)/(-1.427*0.13449)
      } else if(age >=16 & age <16.5){
        output <- ((bmi/20.155)^(-1.434)-1)/(-1.434*0.13299)
      } else if(age >=16.5 & age <17){
        output <- ((bmi/20.347)^(-1.438)-1)/(-1.438*0.13177)
      } else if(age >=17 & age <17.5){
        output <- ((bmi/20.512)^(-1.438)-1)/(-1.438*0.13096)
      } else if(age >=17.5 & age <18){
        output <- ((bmi/20.658)^(-1.432)-1)/(-1.432*0.13050)
      } else if(age >=18 & age){
        output <- ((bmi/20.792)^(-1.423)-1)/(-1.423*0.13033)
      }# End if
    }# End if
  }# End if
  
  return(output)
  
}# End function