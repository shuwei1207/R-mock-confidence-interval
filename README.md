我們要製作從andy模型出來的模擬資料，是否吻合原本資料的marginal effect of advert信賴區間，也就是g(b3,b4)= (1-b3)/2*b4

1. 觀察資料，先找出真實的g(b3,b4)= (1-b3)/2*b4，我們命其為cri
2. 建立1000筆資料的模擬資料，令其為mod3，其1000筆資料的上界跟下界分別存在名為upbg[i]跟lowbg[i]的變數中，因此我們得到1000組95%的信心水準
3. 判斷1000組的信心區間是否有包住真實的cri，若cri在區間內設為1，在區間外設為0，用result印出結果
4. 最終的結果為952，意即1000組的模擬資料中有952組可以成功包住cri，成功機率95.2%，與原先的設定一致!