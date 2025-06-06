// tslint:disable
/**
 * SOFI
 * SOFI Sekvensanalyseplatform
 *
 * The version of the OpenAPI document: 0.1.0
 * 
 *
 * NOTE: This class is auto generated by OpenAPI Generator (https://openapi-generator.tech).
 * https://openapi-generator.tech
 * Do not edit the class manually.
 */

import { exists, mapValues } from '../runtime';
/**
 * 
 * @export
 * @interface Exception
 */
export interface Exception  {
    /**
     * 
     * @type {string}
     * @memberof Exception
     */
    source?: string;
    /**
     * 
     * @type {string}
     * @memberof Exception
     */
    message: string;
    /**
     * 
     * @type {string}
     * @memberof Exception
     */
    stackTrace?: string;
    /**
     * 
     * @type {object}
     * @memberof Exception
     */
    data?: object;
    /**
     * 
     * @type {string}
     * @memberof Exception
     */
    targetSite?: string;
    /**
     * 
     * @type {string}
     * @memberof Exception
     */
    helpLink?: string;
    /**
     * 
     * @type {Exception}
     * @memberof Exception
     */
    innerException?: Exception;
}

export function ExceptionFromJSON(json: any): Exception {
    return {
        'source': !exists(json, 'source') ? undefined : json['source'],
        'message': json['message'],
        'stackTrace': !exists(json, 'stackTrace') ? undefined : json['stackTrace'],
        'data': !exists(json, 'data') ? undefined : json['data'],
        'targetSite': !exists(json, 'targetSite') ? undefined : json['targetSite'],
        'helpLink': !exists(json, 'helpLink') ? undefined : json['helpLink'],
        'innerException': !exists(json, 'innerException') ? undefined : ExceptionFromJSON(json['innerException']),
    };
}

export function ExceptionToJSON(value?: Exception): any {
    if (value === undefined) {
        return undefined;
    }
    return {
        'source': value.source,
        'message': value.message,
        'stackTrace': value.stackTrace,
        'data': value.data,
        'targetSite': value.targetSite,
        'helpLink': value.helpLink,
        'innerException': ExceptionToJSON(value.innerException),
    };
}


