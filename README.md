# Planetary-motion-3.5
## Описание задачи
Проверить что полученная в Planetary-motion-3.1 траектория движения планеты является эллипсом или гиперболой. Получить параметры эллипса. Для этого необходимо реализовать функцию `curve2order`: `[a,b, x0, y0, cos(beta), type] = curve2order(x, y)`. \
Входные данные - массив координат планет. \
Выходные данные - большая `a` и малая `b` полуоси, координаты `(x0, y0)` центра кривой второго порядка, `cos(beta)` косинус угла наклона большой полуоси к оси _x_, и тип кривой (**+1** - эллипс, **-1** - гипербола).
