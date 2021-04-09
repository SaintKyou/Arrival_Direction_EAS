#ifndef STRUCTS_H
#define STRUCTS_H

// EAS parameters: Teta, Fi, XAxisShift, YAxisShift, NeNKGlong, sNKGlong
// NEVOD-EAS data: EdepDetNE[9][4][4], TimDetNE[9][4][4][4]

#pragma pack(push,1)
struct TDateTimeKadr
{
    union{
        struct{
            unsigned char hsecond; //сотые секунды
            unsigned char second;  //секунда
            unsigned char minute;  //минута
            unsigned char hour; //час
        }tm;
        unsigned long time;
    };
    union{
        struct{
            unsigned char day; //день
            unsigned char month; //месяц
            unsigned short year; //год
        }dt;
        unsigned long date;
    };
};
struct HEADER
{
    char start[6]; //  слово "start" - метка начала записи.
    char id[38]; // идентификатор события
    unsigned int lendata; // Длина следующих за заголовком данных в байтах.
};
struct SCOORD
{
    float x;
    float y;
    float z;
};
struct STATION
{
    short maska;//маска отключенных счетчиков ДС
    float Amp, Q; // амплитуда и заряд станции.
    float Amp_add, Q_add; // амплитуда и заряд дополнительного ФЭУ.
    float time, time_add; // ns от начала буфера.
};
struct CLUSTER
{
    char id[38]; // идентификатор события:
//    len(id) = 36: 2018-12-25_1_e_0_0:00:35.375.341.730\x00\x00
//    len(id) = 37: 2018-12-25_1_e_0_10:00:35.375.341.730\x00
    short maska; // маска неработающих станций кластера.
    short hit; // маска сработавших станций.
    long time; // наносекунды.
    STATION st[4]; // массив с данными станций.
};
struct DATA
{
    short maska; // маска неработающих кластеров.
    short Hit; // маска сработавших кластеров.
    CLUSTER clust[9]; // массив с данными кластеров.
};
struct SSETUP
{
    TDateTimeKadr t; // время записи данных.
    short NClust; // число установленных кластеров в установке.
    int maska[4]; // битовая маска отключённых ДС в кластерах (индекс массива соответствует номеру ДС)
    short TRGcase[9]; // триггерное условие (количества совпадений).
    short TRGgate[9]; // ворота для совпадений в кластере 120 нс.
    short TRGtime[9]; // время формирования триггера в буфере 464*5=2320 нс.
    short OpticDelay[9]; // задержка сигнала синхронизации.
    SCOORD point[9][4]; // пространственные координаты детектирующих станций.
};
struct SQPEAK
{
    TDateTimeKadr dt1;//дата и время начала набора данных RUN_begin
    TDateTimeKadr dt2;//дата и время окончания набора данных RUN_end
    short NClust;//	число кластеров в установке
    float Temperature; // средняя за набор температура.
    float SigT; // среднеквадратичное отклонение температуры за набор.
    float Qpic[9][4]; // наиболее вероятная калибровочная амплитуда.
    float SQpic[9][4]; // погрешность наиболее вероятной амплитуды.
    float PeakToVal[9][4]; // отношение пик/долина.
};
struct SCROSSLINK
{
    TDateTimeKadr t1; // время начала данных.
    TDateTimeKadr t2; // время окончания данных.
    short NClust;//	число кластеров в установке
    float CrossLink[9][4]; // коэффициенты сшивки.
    float SCrossLink[9][4]; // погрешность коэффициентов сшивки.
};
struct SMONIT
{
    TDateTimeKadr t1; // время начала интервала.
    TDateTimeKadr t2; // время окончания интервала.
    short NClust;//	число кластеров в установке
    int maska[4]; // маска отключённых ДС в кластерах.
    float Press; // среднее атмосферное давление.
    float Temperature; // средняя температура окружающей среды.
    float TempCl[9]; // средняя температура внутри локальных пунктов кластеров.
    float Rate[9][5]; // средний темп счета кластера и его ДС.
    float SRate[9][5]; // среднеквадратичное отклонение темпа счета кластера и его ДС.
    float Pds[9][4][2]; // среднее значение пьедесталов измерительных каналов кластера.
    float SPds[9][4][2]; // среднеквадратичное отклонение пьедесталов измерительных каналов кластера.
};
struct SDELAYS
{
    TDateTimeKadr dt1;//дата и время начала набора данных RUN_begin
    TDateTimeKadr dt2;//дата и время окончания набора данных RUN_end
    short NClust;//	число кластеров в установке
    float Delay[9][6];
};
#pragma pack(pop)

#endif // STRUCTS_H
